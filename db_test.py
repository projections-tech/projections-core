import logging
import logging.config
import psycopg2
from string import ascii_lowercase

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('db_test')

# This is the set of test commands to test tree compatibility with psycopg2


class Tree():
    def __init__(self, name, type, root=None):
        self.name = name
        self.root = root
        self.type = type
        self.uri = None
        self.children = []
        self.meta_parent_path = None


root = Tree('/', 'dir')

for i in range(3):
    child = Tree(ascii_lowercase[i], 'dir', root=root)
    root.children.append(child)
    for j in range(3):
        sub_child = Tree(name=ascii_lowercase[j], type='dir', root=child)
        child.children.append(sub_child)
        for k in range(2):
            sub_sub_child = Tree(name=ascii_lowercase[k] + '.txt', type='file', root=sub_child)
            sub_child.children.append(sub_sub_child)
            sub_sub_child_meta = Tree(name=ascii_lowercase[k] + '.txt.meta', type='meta', root=sub_child)
            sub_sub_child_meta.meta_parent_path = ['/'] + [ascii_lowercase[el] for el in [i, j, k]]
            sub_sub_child_meta.meta_parent_path[-1] += '.txt'
            sub_child.children.append(sub_sub_child_meta)



def db_build_tree(cursor, tree_table, metadata_table, root, parent_id=None):
    cursor.execute("SELECT relname FROM pg_class WHERE relname = '{}';".format(TREE_TABLE_NAME))
    is_tree_table_exist = cursor.fetchone()
    if not bool(is_tree_table_exist):
        cursor.execute("CREATE TABLE {} (node_id serial PRIMARY KEY, parent_id integer, name varchar,"
                            " path varchar[][], type varchar);".format(tree_table))

        cursor.execute("CREATE TABLE {} (meta_id serial PRIMARY KEY, node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE,"
                    " parent_node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE,"
                    " meta_parent_path varchar);".format(metadata_table))

    name = root.name
    projection_type = root.type

    insertion_command = "INSERT INTO {0} (parent_id, name, type) VALUES (%s, %s, %s) RETURNING node_id".format(tree_table)

    cursor.execute(insertion_command, (parent_id, name, projection_type))
    parent_id = cursor.fetchone()

    if root.type == 'meta':
        insertion_command = 'INSERT INTO {0} (node_id, meta_parent_path) VALUES (%s, %s)'.format(metadata_table)
        cursor.execute(insertion_command, (parent_id, root.meta_parent_path))

    for child in root.children:
        db_build_tree(cursor, tree_table, metadata_table, child, parent_id)


def db_fill_metadata_table(cursor, table_1, table_2):

    insertion_command = """
    UPDATE {1} SET parent_node_id={0}.node_id FROM {0} WHERE {0}.path=meta_parent_path::varchar[];
    """.format(table_1, table_2)

    cursor.execute(insertion_command)


def db_update_tree_paths(cursor, table):

    recursive_update_command = """

    WITH RECURSIVE ancestors_table AS (
        SELECT node_id, ARRAY[]::integer[] AS ancestors, ARRAY[{0}.name]::varchar[] AS anc_paths
        FROM {0} WHERE parent_id IS NULL

        UNION ALL

        SELECT {0}.node_id, ancestors_table.ancestors || {0}.parent_id, ancestors_table.anc_paths || {0}.name
        FROM {0}, ancestors_table
        WHERE {0}.parent_id = ancestors_table.node_id
    )
    UPDATE {0} SET path=ancestors_table.anc_paths FROM ancestors_table WHERE {0}.node_id=ancestors_table.node_id;

    """.format(table)

    cursor.execute(recursive_update_command)


def db_update_node_descendants_paths(cursor, tree_table, metadata_table, node_id, new_parent_node_id):

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

    """.format(tree_table, new_parent_node_id)

    cursor.execute(recursive_update_command)

    metadata_path_update_command = """
    UPDATE {1} SET meta_parent_path={0}.path FROM {0} WHERE parent_node_id = {2} AND {0}.node_id = {2}
    """.format(tree_table, metadata_table, node_id)
    logger.debug('Command : {}'.format(metadata_path_update_command))
    cursor.execute(metadata_path_update_command)

def db_move_node(cursor, tree_table, metadata_table, node_path, new_node_root_path):
    assert isinstance(node_path, list), 'Input old path is not a list!'
    assert isinstance(new_node_root_path, list), 'Input new path is not a list!'

    node_ids_command = """
    WITH
    old_root_node AS (
        SELECT node_id AS old_root_node_id FROM {0} WHERE path=%s::varchar[]
    ),
    new_root_node AS (
        SELECT node_id AS new_root_node_id FROM {0} WHERE path=%s::varchar[]
    )
    SELECT * FROM old_root_node, new_root_node;

    """.format(tree_table)

    cursor.execute(node_ids_command, (node_path, new_node_root_path))

    fetch_result = cursor.fetchone()

    if fetch_result is not None:
        node_id, new_parent_id = fetch_result
    else:
        raise OSError('Attempting to move non existant node on path: {}'.format(node_path))

    move_command = """
    UPDATE {0} SET parent_id={1} WHERE node_id = {2}
    """.format(tree_table, new_parent_id, node_id)
    cursor.execute(move_command)

    db_update_node_descendants_paths(cursor, tree_table, metadata_table, node_id, new_parent_id)


def db_remove_node(cursor, table, node_path):
    """
    :param cursor:
    :param table:
    :param node_id:
    :return:
    """
    assert isinstance(node_path, list), 'Node path is not a list!'

    fetch_node_to_remove_id_command = """
    SELECT node_id FROM {0} WHERE path=%s::varchar[]
    """.format(table)
    cursor.execute(fetch_node_to_remove_id_command, (node_path, ) )

    fetch_result = cursor.fetchone()

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
    """.format(table, node_to_remove_id)

    cursor.execute(node_removal_command)

TREE_TABLE_NAME = 'tree_table'
META_TABLE_NAME = 'metadata_table'

with psycopg2.connect("dbname=test user=viktor") as conn:
    with conn.cursor() as cursor:

        cursor.execute("SELECT relname FROM pg_class WHERE relname = '{}';".format(TREE_TABLE_NAME))
        is_tree_table_exist = cursor.fetchone()
        if not bool(is_tree_table_exist):
            db_build_tree(cursor, TREE_TABLE_NAME, META_TABLE_NAME, root)
            db_update_tree_paths(cursor, TREE_TABLE_NAME)
            db_fill_metadata_table(cursor, TREE_TABLE_NAME, META_TABLE_NAME)

        db_move_node(cursor=cursor, tree_table=TREE_TABLE_NAME, metadata_table=META_TABLE_NAME, node_path=['/', 'a', 'a', 'b.txt'], new_node_root_path=['/'])

        cursor.execute('SELECT * FROM {}'.format(TREE_TABLE_NAME))
        for el in cursor:
            logger.debug(el)



        cursor.execute('SELECT * FROM {}'.format(TREE_TABLE_NAME))
        for el in cursor:
            logger.debug(el)

        conn.commit()

        cursor.close()

conn.close()