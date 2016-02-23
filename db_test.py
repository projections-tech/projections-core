import psycopg2
from string import ascii_lowercase

# This is the set of test commands to test tree compatibility with psycopg2


class Tree():
    def __init__(self, name, root=None):
        self.name = name
        self.root = root
        self.children = []


root = Tree('/')

for i in range(3):
    child = Tree(ascii_lowercase[i], root=root)
    root.children.append(child)
    for j in range(3):
        sub_child = Tree(name=ascii_lowercase[j], root=child)
        child.children.append(sub_child)
        for k in range(2):
            sub_sub_child = Tree(name=ascii_lowercase[k] + '.txt', root=sub_child)
            sub_child.children.append(sub_sub_child)


def db_build_tree(cursor, table, root, parent_id=None):
    name = root.name

    cur.execute("INSERT INTO test (parent_id, name) VALUES (%s, %s) RETURNING node_id", (parent_id, name))
    parent_id = cur.fetchone()
    for child in root.children:
        db_build_tree(cursor, table, child, parent_id)


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

    cur.execute(recursive_update_command)


def db_update_node_descendants_paths(cursor, table, node_id):

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

    """.format(table, node_id)

    cur.execute(recursive_update_command)


def db_move_node(cursor, table, node_path, new_node_root_path):
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

    """.format(table)

    cursor.execute(node_ids_command, (node_path, new_node_root_path))
    node_id, new_parent_id = cursor.fetchone()

    move_command = """
    UPDATE {0} SET parent_id={1} WHERE node_id = {2}
    """.format(table, new_parent_id, node_id)
    cursor.execute(move_command)

    db_update_node_descendants_paths(cursor, table, new_parent_id)

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
        raise OSError('Attempting to delete non existant node!')

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


with psycopg2.connect("dbname=test user=enmce") as conn:
    with conn.cursor() as cur:

        cur.execute("SELECT relname FROM pg_class WHERE relname = 'test';")
        is_table_exist = cur.fetchone()

        if not bool(is_table_exist):
            cur.execute("CREATE TABLE test (node_id serial PRIMARY KEY, parent_id integer, name varchar, path varchar[][]);")
            db_build_tree(cur, 'test', root)
            db_update_tree_paths(cur, 'test')

        cur.execute('SELECT * FROM test')
        for el in cur:
            print(el)

        db_move_node(cur, 'test', ['/', 'a'], ['/', 'b'])

        cur.execute('SELECT * FROM test')
        for el in cur:
            print(el)


        conn.commit()

        cur.close()

conn.close()