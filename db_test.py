import logging
import logging.config
import psycopg2
import psycopg2.extras
from string import ascii_lowercase

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('db_test')


class Tree:
    def __init__(self, name, type, root=None):
        self.name = name
        self.root = root
        self.type = type
        self.uri = None
        self.children = {}
        self.meta_parent_path = None
        self.meta_contents = None


root = Tree('/', 'dir')

for i in range(3):
    child = Tree(ascii_lowercase[i], 'dir', root=root)
    root.children[i] = child
    for j in range(3):
        sub_child = Tree(name=ascii_lowercase[j], type='dir', root=child)
        child.children[j] = sub_child
        for k in range(2):
            sub_sub_child = Tree(name=ascii_lowercase[k] + '.txt', type='file', root=sub_child)
            sub_child.children[k] = sub_sub_child
            sub_sub_child_meta = Tree(name=ascii_lowercase[k] + '.txt.meta', type='meta', root=sub_child)
            sub_sub_child_meta.meta_parent_path = ['/'] + [ascii_lowercase[el] for el in [i, j, k]]
            sub_sub_child_meta.meta_parent_path[-1] += '.txt'
            sub_sub_child_meta.meta_contents = {'valid': True, 'name': str(k)}
            if k == 1:
                sub_sub_child_meta.meta_contents['valid'] = False
            sub_child.children[str(k)+'meta'] = sub_sub_child_meta


class DBProjector:
    """
    This class contains test interface to PostgresSQL database, which will be used to store and work with projections.
    Current implementation is subject of future rewrites and optimisations,
    """

    def __init__(self, user_name, database_name, prototypes_tree):
        """
        This method initializes database which will hold projection tree and associated metadata
        """
        self.tree_table_name = 'tree_table'
        self.metadata_table_name = 'metadata_table'

        # Initializing psycopg JSONB support
        psycopg2.extras.register_json(oid=3802, array_oid=3807)
        # Opening connection with database
        self.db_connection = psycopg2.connect("dbname={db_name} user={user_name}".format(db_name=database_name,
                                                                                         user_name=user_name))
        # Creating cursor, which will be used to interact with database
        self.cursor = self.db_connection.cursor()

        # Creating tables with which DBProjector will work
        self.db_create_tables()
        # Building projections tree
        self.db_build_tree(prototypes_tree)
        # Updating paths of nodes in tree_table

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
                                "path varchar[][], "
                                "type varchar"
                                ");".format(self.tree_table_name))

            self.cursor.execute("CREATE TABLE {} ("
                                "meta_id serial PRIMARY KEY,"
                                "node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE,"
                                "parent_node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE,"
                                "meta_parent_path varchar, "
                                "meta_contents jsonb);".format(self.metadata_table_name))

    def db_build_tree(self, prototype_tree, parent_id=None, projection_path=None):
        """
        This method recursively builds projections file structure in tree_table from prototype tree
        :param prototype_tree: Prototype Tree object
        :param parent_id: node_id of parent projection in tree_table
        """
        #TODO add actual driver support
        #TODO solve Metadata-Projection binding!

        projection_name = prototype_tree.name
        projection_type = prototype_tree.type
        if not projection_path is None:
            # Copying projection path before current path appending
            projection_path = projection_path[:]
            projection_path.append(projection_name)
        else:
            projection_path = [projection_name]

        insertion_command = "INSERT INTO {0} (parent_id, name, type, path) " \
                            "SELECT %(p_id)s, %(name)s, %(type)s, %(path)s" \
                            "WHERE NOT EXISTS (" \
                            "SELECT * FROM {0} " \
                            "WHERE path=%(path)s::varchar[])" \
                            "RETURNING node_id".format(self.tree_table_name)

        self.cursor.execute(insertion_command, {'p_id': parent_id,
                                                'name': projection_name,
                                                'type': projection_type,
                                                'path': projection_path})
        parent_id = self.cursor.fetchone()

        if prototype_tree.type == 'meta':
            insertion_command = 'INSERT INTO {0} (node_id, meta_parent_path, meta_contents) ' \
                                'VALUES (%s, %s, %s)'.format(self.metadata_table_name)
            self.cursor.execute(insertion_command, (parent_id, prototype_tree.meta_parent_path,
                                                    psycopg2.extras.Json(prototype_tree.meta_contents)))

        for _, child in prototype_tree.children.items():
            self.db_build_tree(child, parent_id, projection_path)

        self.db_connection.commit()

    def db_fill_metadata_table(self):
        """
        This method does object-metadata binding for all objects in tree_table which have metadata
        """
        # Forming insertion command
        insertion_command = """
        UPDATE {1} SET parent_node_id={0}.node_id FROM {0} WHERE {0}.path=meta_parent_path::varchar[];
        """.format(self.tree_table_name, self.metadata_table_name)
        # Executing and committing result on success
        self.cursor.execute(insertion_command)
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

        metadata_path_update_command = """
        WITH RECURSIVE

        descendants_table AS (
            WITH RECURSIVE tree AS (
                SELECT node_id, ARRAY[{2}]::integer[] AS ancestors
                FROM {0} WHERE parent_id = {2}

                UNION ALL

                SELECT {0}.node_id, tree.ancestors || {0}.parent_id
                FROM {0}, tree
                WHERE {0}.parent_id = tree.node_id
            )
            SELECT node_id FROM tree WHERE {2} = ANY(tree.ancestors)
        ),

        ancestors_table AS (
            SELECT node_id, ARRAY[]::integer[] AS ancestors, ARRAY[{0}.name]::varchar[] AS anc_paths
            FROM {0} WHERE parent_id IS NULL

            UNION ALL

            SELECT {0}.node_id, ancestors_table.ancestors || {0}.parent_id, ancestors_table.anc_paths || {0}.name
            FROM {0}, ancestors_table
            WHERE {0}.parent_id = ancestors_table.node_id
        )

        UPDATE {1} SET meta_parent_path={0}.path
        FROM {0}, descendants_table
        WHERE parent_node_id = descendants_table.node_id AND {0}.node_id = descendants_table.node_id
        """.format(self.tree_table_name, self.metadata_table_name, moved_node_id)

        self.cursor.execute(metadata_path_update_command)
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
            SELECT node_id, path AS target_id FROM {0} WHERE path=%s::varchar[]
        ),
        metadata AS (
            SELECT node_id AS metadata_id FROM {0} WHERE path=%s::varchar[]
        )
        INSERT INTO {1} (parent_node_id, meta_parent_path, node_id) SELECT * FROM target, metadata RETURNING node_id
        """.format(self.tree_table_name, self.metadata_table_name)
        self.cursor.execute(binding_command, (target_path, metadata_path))
        self.db_connection.commit()


# Smoke testing commences here
if __name__ == '__main__':
    test_db_projector = DBProjector('viktor', 'test', root)


