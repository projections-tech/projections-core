import copy
import logging
import logging.config
import os
import types

import objectpath
import psycopg2
import psycopg2.extras

import iontorrent
from projections import PrototypeDeserializer
from tests.mock import MockResource

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('db_test')


class DBProjector:
    """
    This class contains test interface to PostgresSQL database, which will be used to store and work with projections.
    Current implementation is subject of future rewrites and optimisations,
    """

    def __init__(self, projection_driver, user_name, database_name, prototypes_tree):
        """
        This method initializes database which will hold projection tree and associated metadata
        """
        self.tree_table_name = 'tree_table'
        self.metadata_table_name = 'metadata_table'
        self.prototypes_tree = prototypes_tree

        # Initializing projection driver
        self.projection_driver = projection_driver

        # Initializing psycopg JSONB support
        psycopg2.extras.register_json(oid=3802, array_oid=3807)
        # Opening connection with database
        self.db_connection = psycopg2.connect("dbname={db_name} user={user_name}".format(db_name=database_name,
                                                                                         user_name=user_name))
        # Creating cursor, which will be used to interact with database
        self.cursor = self.db_connection.cursor()

        # Creating tables with which DBProjector will work
        self.db_create_tables()

    def __del__(self):
        """
        This method closes cursor and database connection at call to destructor method
        """
        self.cursor.close()
        self.db_connection.close()

    def db_create_tables(self):
        """
        This method creates tables which will be used to store projections structure and associated metadata
        """
        # Checking if tree_table exists in pg_class table
        self.cursor.execute("SELECT relname FROM pg_class WHERE relname = '{}';".format(self.tree_table_name))
        is_tree_table_exist = self.cursor.fetchone()

        # If no tree table found, create required tables
        if not bool(is_tree_table_exist):
            self.cursor.execute("CREATE TABLE {} ("
                                "node_id serial PRIMARY KEY, "
                                "parent_id integer, "
                                "name varchar, "
                                "uri varchar, "
                                "path varchar[][], "
                                "type varchar"
                                ");".format(self.tree_table_name))

            self.cursor.execute("DROP TABLE IF EXISTS {0}".format(self.metadata_table_name))

            self.cursor.execute("CREATE TABLE {} ("
                                "meta_id serial PRIMARY KEY,"
                                "node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE,"
                                "parent_node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE,"
                                "meta_contents jsonb);".format(self.metadata_table_name))

            self.db_build_tree({'/': self.prototypes_tree}, projection_path=['/'],
                           root_projection_uri='experiment?status=run&limit=1&order_by=-id')

    def db_build_tree(self, prototypes, parent_id=None, projection_path=None, root_projection_uri=None,
                      parent_id_for_meta=None):
        """
        This method recursively builds projections file structure in tree_table from prototype tree
        :param prototypes: Prototype Tree object
        :param parent_id: node_id of parent projection in tree_table
        :param projection_path:
        :param root_projection_uri
        """
        #TODO solve Metadata-Projection binding!

        environment = None
        content = None
        path = os.path

        # This is environment in which projections are created (parent_projection content)
        # TODO: in many cases it means double request to parent projection resource so it should be optimized
        # We don`t want to change driver contents, hence we made deep copy of dict
        environment = copy.deepcopy(self.projection_driver.get_uri_contents_as_dict(root_projection_uri))

        # For every prototype in collection try to create corresponding projections

        for key, prototype in prototypes.items():
            # Set current prototype context to current environment for children node to use
            prototype.context = environment
            # Get context of current node from contexts of parent nodes
            context = prototype.get_context()
            context = context[::-1]

            # Adding context of upper level prototypes for lower level projections to use
            environment['context'] = context

            # Creating tree of environment contents which will be parsed by ObjectPath
            tree = objectpath.Tree(environment)
            URIs = tree.execute(prototype.uri)
            # Object path sometimes returns generator if user uses selectors, for consistency expand it using
            # list comprehension
            if isinstance(URIs, types.GeneratorType):
                URIs = [el for el in URIs]
            # Treating URIs as list for consistency
            if not isinstance(URIs, list):
                URIs = [URIs]

            # We get projection URIs based on environment and prototype properties
            # Every URI corresponds to projection object
            for uri in URIs:
                current_projection_path = projection_path[:]
                # Get content for a projection
                # We don`t want to change driver contents, hence we made deep copy
                content = copy.deepcopy(self.projection_driver.get_uri_contents_as_dict(uri))

                # Adding environment to use by prototype
                content['environment'] = environment
                content['context'] = context

                # Creating tree which will be parsed by ObjectPath
                tree = objectpath.Tree(content)
                name = tree.execute(prototype.name)

                # Object path sometimes returns generator if user uses selectors, for consistency expand it using
                # list comprehension
                if isinstance(name, types.GeneratorType):
                    name = [el for el in name]

                if prototype.type != 'metadata':
                    set_parent_id = parent_id
                    current_projection_path.append(name)
                else:
                    current_projection_path = current_projection_path[:-1]
                    current_projection_path.append(name)
                    set_parent_id = prototype.meta_parent_id

                insertion_command = "INSERT INTO {0} (parent_id, name, uri, type, path) " \
                                    "SELECT %(p_id)s, %(name)s, %(uri)s, %(type)s, %(path)s" \
                                    "WHERE NOT EXISTS (" \
                                    "SELECT * FROM {0} " \
                                    "WHERE path=%(path)s::varchar[])" \
                                    "RETURNING node_id".format(self.tree_table_name)

                self.cursor.execute(insertion_command, {'p_id': set_parent_id,
                                                        'name': name,
                                                        'type': prototype.type,
                                                        'path': current_projection_path,
                                                        'uri': uri})
                new_parent_id = self.cursor.fetchone()

                has_meta = False
                for key, child in prototype.children.items():
                    if child.type == 'metadata':
                        child.meta_parent_id = set_parent_id
                        has_meta = True

                if prototype.type == 'metadata' and not (new_parent_id is None and parent_id_for_meta is None):
                    meta_table_insertion_command = "INSERT INTO {0} (node_id, parent_node_id, meta_contents) " \
                                                   "VALUES (%(node_id)s, %(parent_node_id)s, %(meta_contents)s)".format(self.metadata_table_name)

                    metadata_contents = self.projection_driver.get_uri_contents_as_dict(uri)

                    self.cursor.execute(meta_table_insertion_command, {'node_id': new_parent_id,
                                                                       'parent_node_id': parent_id_for_meta,
                                                                       'meta_contents': psycopg2.extras.Json(metadata_contents)})

                if has_meta:
                    self.db_build_tree(prototype.children, new_parent_id, current_projection_path,
                                       root_projection_uri=uri, parent_id_for_meta=new_parent_id)
                else:
                    self.db_build_tree(prototype.children, new_parent_id, current_projection_path,
                                       root_projection_uri=uri)

                self.db_connection.commit()

    def db_update_node_descendants_paths(self, moved_node_id, new_parent_node_id):
        """
        This method updates node descendants paths, and path to parent node in metadata table
        :param moved_node_id: node_id of moved node in tree_table
        :param new_parent_node_id: node_id of new parent node in tree_table
        """

        recursive_update_command = """

        WITH RECURSIVE

        descendants_table AS (
            WITH RECURSIVE tree AS (
                SELECT node_id, ARRAY[{1}]::integer[] AS ancestors
                FROM {0} WHERE parent_id = {1}

                UNION ALL

                SELECT {0}.node_id, tree.ancestors || {0}.parent_id
                FROM {0}, tree
                WHERE {0}.parent_id = tree.node_id
            )
            SELECT node_id FROM tree WHERE {1} = ANY(tree.ancestors)
        ),

        ancestors_table AS (
            SELECT node_id, ARRAY[]::integer[] AS ancestors, ARRAY[{0}.name]::varchar[] AS anc_paths
            FROM {0} WHERE parent_id IS NULL

            UNION ALL

            SELECT {0}.node_id, ancestors_table.ancestors || {0}.parent_id, ancestors_table.anc_paths || {0}.name
            FROM {0}, ancestors_table
            WHERE {0}.parent_id = ancestors_table.node_id
        )

        UPDATE {0} SET path=ancestors_table.anc_paths FROM ancestors_table, descendants_table
        WHERE {0}.node_id=ancestors_table.node_id AND {0}.node_id=descendants_table.node_id;

        """.format(self.tree_table_name, new_parent_node_id)

        self.cursor.execute(recursive_update_command)
        self.db_connection.commit()

    def db_move_node(self, node_to_move_path, new_node_root_path):
        """
        This methods moves node in a tree to a new parent node, and updates node paths accordingly
        :param node_to_move_path: path to node which will be moved as list of strings
        :param new_node_root_path: path to new parent node list of strings
        """
        assert isinstance(node_to_move_path, list), 'Input old path is not a list!'
        assert isinstance(new_node_root_path, list), 'Input new path is not a list!'
        if node_to_move_path == new_node_root_path:
            raise Exception('Cannot move node onto itself!')

        node_ids_command = """
        WITH
        old_root_node AS (
            SELECT node_id AS old_root_node_id FROM {0} WHERE path=%s::varchar[]
        ),
        new_root_node AS (
            SELECT node_id AS new_root_node_id FROM {0} WHERE path=%s::varchar[]
        )
        SELECT * FROM old_root_node, new_root_node;

        """.format(self.tree_table_name)

        self.cursor.execute(node_ids_command, (node_to_move_path, new_node_root_path))

        fetch_result = self.cursor.fetchone()

        if fetch_result is not None:
            node_id, new_parent_id = fetch_result
        else:
            raise OSError('Attempting to move non existant node on path: {}'.format(node_to_move_path))

        move_command = """
        UPDATE {0} SET parent_id={1} WHERE node_id = {2}
        """.format(self.tree_table_name, new_parent_id, node_id)
        self.cursor.execute(move_command)
        self.db_connection.commit()

        self.db_update_node_descendants_paths(node_id, new_parent_id)

    def db_remove_node(self, node_path):
        """
        This method removes node from tree_table and all of it`s descendants
        :param node_path: path to node list of strings
        """
        assert isinstance(node_path, list), 'Node path is not a list!'

        fetch_node_to_remove_id_command = """
        SELECT node_id FROM {0} WHERE path=%s::varchar[]
        """.format(self.tree_table_name)
        self.cursor.execute(fetch_node_to_remove_id_command, (node_path, ))
        self.db_connection.commit()

        fetch_result = self.cursor.fetchone()

        if fetch_result is not None:
            node_to_remove_id = fetch_result[0]
        else:
            raise OSError('Attempting to delete non existant node on path {}'.format(node_path))

        node_removal_command = """
        WITH RECURSIVE descendants_table AS (
            SELECT node_id, ARRAY[]::integer[] AS ancestors
            FROM {0} WHERE parent_id IS NULL

            UNION ALL

            SELECT {0}.node_id, descendants_table.ancestors || {0}.parent_id
            FROM {0}, descendants_table
            WHERE {0}.parent_id = descendants_table.node_id
        )
        DELETE FROM {0}
        WHERE node_id IN (
            SELECT node_id FROM descendants_table WHERE {1} = ANY(descendants_table.ancestors)
        ) OR node_id = {1};
        """.format(self.tree_table_name, node_to_remove_id)

        self.cursor.execute(node_removal_command)
        self.db_connection.commit()

    def db_bind_metadata_to_path(self, target_path, metadata_path):
        """
        This method binds metadata object on path to target path
        :param target_path: path to target with which metadata is binded, list of stings
        :param metadata_path: path to metadata object which will be binded, list of strings
        """

        binding_command = """
        WITH
        target AS (
            SELECT node_id AS target_id FROM {0} WHERE path=%s::varchar[]
        ),
        metadata AS (
            SELECT node_id AS metadata_id FROM {0} WHERE path=%s::varchar[]
        )
        INSERT INTO {1} (parent_node_id, node_id) SELECT target.target_id, metadata.metadata_id FROM target, metadata
        WHERE NOT EXISTS (SELECT * FROM target, metadata, metadata_table
        WHERE metadata_table.node_id = metadata.metadata_id AND metadata_table.parent_node_id = target.target_id)
        """.format(self.tree_table_name, self.metadata_table_name)
        self.cursor.execute(binding_command, (target_path, metadata_path))

        self.db_connection.commit()

# Smoke testing commences here
if __name__ == '__main__':
    USER = 'user'
    PASSWORD = 'password'
    mock_resource = MockResource('tests/torrent_suite_mock.json')
    projection_configuration = PrototypeDeserializer('test_full_torrent_suite_config.yaml')
    driver = iontorrent.TorrentSuiteDriver(projection_configuration.resource_uri, USER, PASSWORD)

    test_db_projector = DBProjector(driver, 'viktor', 'test',
                                    projection_configuration.prototype_tree)

    test_db_projector.db_move_node(['/', 'test_experiment_1_plannedexperiment.meta'], ['/', 'test_experiment_2'])

    mock_resource.deactivate()


