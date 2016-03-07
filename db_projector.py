import copy
import getpass
import io
import logging
import logging.config
import os
import stat
import time
import types

import objectpath
import psycopg2
import psycopg2.extras

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('db_test')


class DBProjector:
    """
    This class contains test interface to PostgresSQL database, which will be used to store and work with projections.
    Current implementation is subject of future rewrites and optimisations,
    """

    def __init__(self, projection_name, projection_driver, prototypes_tree, root_uri):
        """
        This method initializes database which will hold projection tree and associated metadata
        """
        self.projection_name = projection_name

        self.tree_table_name = 'tree_table'
        self.metadata_table_name = 'metadata_table'
        self.projections_attributes_table_name = 'projections_attributes_table'
        self.prototypes_tree = prototypes_tree
        self.root_uri = root_uri

        # Initializing projection driver
        self.projection_driver = projection_driver

        # Initializing psycopg JSONB support
        psycopg2.extras.register_json(oid=3802, array_oid=3807)

        try:
            # Opening connection with database
            self.db_connection = psycopg2.connect(
                "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))
        except Exception as e:
            logger.debug(e)
            logger.debug('Database "projections_database" was not found! Creating database.')
            with psycopg2.connect(dbname='postgres', user=getpass.getuser()) as temp_connection:
                from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
                temp_connection.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)

                with temp_connection.cursor() as temp_cursor:
                    temp_cursor.execute("CREATE DATABASE projections_database")

        # Creating cursor, which will be used to interact with database
        self.cursor = self.db_connection.cursor()

        # Creating tables with which DBProjector will work
        self.db_create_tables()
        logger.debug('Initialized projection: {}'.format(self.projection_name))

    def __del__(self):
        """
        This method closes cursor and database connection at call to destructor method
        """
        self.cursor.close()
        self.db_connection.close()

    def __split_path(self, path):
        """
        Split path string in parts preserving path root
        :param path: path string
        :return: list of path parts strings
        """
        temp = []
        while True:
            head, tail = os.path.split(path)
            if head == '/' or head == '':
                if tail == '':
                    return [head]
                else:
                    temp.append(tail)
                    if not head == '':
                        temp.append(head)
                    return temp[::-1]
            temp.append(tail)
            path = head

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
                                "projection_name varchar,"
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

            self.cursor.execute("DROP TABLE IF EXISTS {0}".format(self.projections_attributes_table_name))

            self.cursor.execute("CREATE TABLE {} ("
                                "attr_id serial PRIMARY KEY,"
                                "node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE, "
                                "st_atime int, "
                                "st_mtime int, "
                                "st_ctime int, "
                                "st_size int, "
                                "st_mode varchar, "
                                "st_nlink int, "
                                "st_ino int);".format(self.projections_attributes_table_name))

        # Add root node to tables
        tree_table_insertion_command = "INSERT INTO {0} (projection_name, parent_id, name, uri, type, path) " \
                                       "SELECT %(proj_name)s, %(p_id)s, %(name)s, %(uri)s, %(type)s, %(path)s " \
                                       "WHERE NOT EXISTS (" \
                                       "SELECT * FROM {0} " \
                                       "WHERE (path=%(path)s::varchar[]) AND projection_name=%(proj_name)s)" \
                                       "RETURNING node_id".format(self.tree_table_name)

        self.cursor.execute(tree_table_insertion_command, {'proj_name': self.projection_name,
                                                           'p_id': None,
                                                           'name': '/',
                                                           'type': 'directory',
                                                           'path': ['/'],
                                                           'uri': self.root_uri})

        new_parent_id = self.cursor.fetchone()

        now = time.time()
        projection_attributes_insertion_command = """
        INSERT INTO {0} (node_id, st_atime, st_mtime, st_ctime, st_size, st_mode, st_nlink, st_ino)
        VALUES (%(node_id)s, %(time)s, %(time)s, %(time)s, %(size)s, %(mode)s, %(nlink)s, %(ino)s);
        """.format(self.projections_attributes_table_name)
        self.cursor.execute(projection_attributes_insertion_command,
                            {'node_id': new_parent_id,
                             'time': now,
                             'size': 0,
                             'mode': 'directory',
                             'nlink': 0,
                             'ino': 1
                             })

        self.db_build_tree({'/': self.prototypes_tree}, projection_path=['/'],
                           root_projection_uri=self.root_uri, parent_id=new_parent_id)

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

                tree_table_insertion_command = "INSERT INTO {0} (projection_name, parent_id, name, uri, type, path) " \
                                               "SELECT %(proj_name)s, %(p_id)s, %(name)s, " \
                                               "%(uri)s, %(type)s, %(path)s" \
                                               "WHERE NOT EXISTS (" \
                                               "SELECT * FROM {0} " \
                                               "WHERE (path=%(path)s::varchar[] AND projection_name=%(proj_name)s))" \
                                               "RETURNING node_id".format(self.tree_table_name)

                self.cursor.execute(tree_table_insertion_command, {'proj_name': self.projection_name,
                                                                   'p_id': set_parent_id,
                                                                   'name': name,
                                                                   'type': prototype.type,
                                                                   'path': current_projection_path,
                                                                   'uri': uri})
                new_parent_id = self.cursor.fetchone()

                if new_parent_id is not None:
                    now = time.time()

                    projection_attributes_insertion_command = """
                    INSERT INTO {0} (node_id, st_atime, st_mtime, st_ctime, st_size, st_mode, st_nlink, st_ino)
                    VALUES (%(node_id)s, %(time)s, %(time)s, %(time)s, %(size)s, %(mode)s, %(nlink)s, %(ino)s);
                    """.format(self.projections_attributes_table_name)
                    self.cursor.execute(projection_attributes_insertion_command,
                                        {'node_id': new_parent_id,
                                         'time': now,
                                         'size': 0,
                                         'mode': prototype.type,
                                         'nlink': 0,
                                         'ino': 1
                                         })

                    has_meta = False
                    for key, child in prototype.children.items():
                        if child.type == 'metadata':
                            child.meta_parent_id = set_parent_id
                            has_meta = True

                    if prototype.type == 'metadata':
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

    def db_update_node_descendants_paths(self, new_parent_node_id):
        """
        This method updates node descendants paths
        :param new_parent_node_id: node_id of new parent node in tree_table which descendants we update
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

    def get_projections_on_path(self, path):
        """
        Method returns list of nodes on path
        :param path: path to node string
        :return: list of nodes paths strings
        """
        path = self.__split_path(path)

        assert isinstance(path, list), 'Path is no a list!'
        paths_request_command = """
        WITH node_on_path AS (
            SELECT node_id AS node_on_path_id FROM {0} WHERE path=%s::varchar[]
        )
        SELECT {0}.name FROM {0}, node_on_path WHERE {0}.parent_id=node_on_path.node_on_path_id
        """.format(self.tree_table_name)

        self.cursor.execute(paths_request_command, (path, ) )
        return [row[0] for row in self.cursor]

    def db_move_projection(self, node_to_move_path, new_node_root_path):
        """
        This methods moves node in a tree to a new parent node, and updates node paths accordingly
        :param node_to_move_path: path to node which will be moved as list of strings
        :param new_node_root_path: path to new parent node list of strings
        """
        node_to_move_path = self.__split_path(node_to_move_path)
        new_node_root_path = self.__split_path(new_node_root_path)

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

        self.db_update_node_descendants_paths(new_parent_id)

    def db_remove_projection(self, node_path):
        """
        This method removes node from tree_table and all of it`s descendants
        :param node_path: path to node list of strings
        """
        node_path = self.__split_path(node_path)

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

        target_path = self.__split_path(target_path)
        metadata_path = self.__split_path(metadata_path)

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

    def is_managing_path(self, path):
        """
        Checj if database is managing path
        :param path: projection path string
        :return: bool
        """
        path = self.__split_path(path)

        # This command checks existance of projection row in tree_table by path
        projection_on_path_existance_check = """
        SELECT EXISTS(
            SELECT * FROM {0} WHERE path=%s::varchar[]
        )
        """.format(self.tree_table_name)
        # Run command and return check result as bool
        self.cursor.execute(projection_on_path_existance_check, (path, ))
        return self.cursor.fetchone()[0]

    def get_attributes(self, path):
        """
        Get attributes of projection on given path
        """

        path = self.__split_path(path)

        attributes_order = ['st_atime', 'st_mtime', 'st_ctime', 'st_size', 'st_mode', 'st_nlink', 'st_ino']

        get_attributes_command = """
        WITH projection_on_path AS (
            SELECT node_id FROM {0} WHERE path = %s::varchar[]
        )
        SELECT {2} FROM {1} JOIN projection_on_path ON {1}.node_id=projection_on_path.node_id
        """.format(self.tree_table_name, self.projections_attributes_table_name, ', '.join(attributes_order))

        self.cursor.execute(get_attributes_command, (path, ))
        # Fetching results of query
        attributes = self.cursor.fetchone()

        # Setting projection attributes dictionary using query results
        attributes = {el[0]: el[1] for el in zip(attributes_order, attributes)}

        # Setting appropriate types access modes for projection types for FUSE
        access_modes = {'file': (stat.S_IFREG | 0o0777),
                       'metadata': (stat.S_IFREG | 0o0777),
                       'directory': (stat.S_IFDIR | 0o0777)}

        attributes['st_mode'] = access_modes[attributes['st_mode']]
        return attributes

    def open_resource(self, path):
        """
        Opens resource on path and returns it`s header and contents stream
        :param path path string
        :return file_header
        :return resource_io
        """
        path = self.__split_path(path)

        projection_id_and_uri_query = """
        SELECT node_id, uri, type FROM {0} WHERE path = %s::varchar[]
        """.format(self.tree_table_name)

        self.cursor.execute(projection_id_and_uri_query, (path, ))

        node_id, uri, projecton_type = self.cursor.fetchone()

        content = self.projection_driver.get_uri_contents_as_bytes(uri)
        logger.info('Got path content: %s\n', path)

        projection_size = len(content)

        update_projection_attributes_command = """
        UPDATE {0} SET st_size=%s WHERE node_id=%s
        """.format(self.projections_attributes_table_name)
        self.cursor.execute(update_projection_attributes_command, (projection_size, node_id,))

        file_header = 3
        resource_io = io.BytesIO(content)

        self.db_connection.commit()

        return file_header, resource_io