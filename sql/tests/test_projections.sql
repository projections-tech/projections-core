-- Create pgTAP extension for testing
CREATE EXTENSION IF NOT EXISTS pgtap;

-- Setup function that can be used to initialize table records before running tests
CREATE OR REPLACE FUNCTION setup_projections() RETURNS SETOF text AS 
$$
BEGIN
	RETURN;
END;
$$ LANGUAGE plpgsql;

-- Setup function that can be used to initialize table records before running tests
CREATE OR REPLACE FUNCTION test_tree_functions() RETURNS SETOF text AS 
$BODY$
DECLARE
    _parent_id bigint;
    _child_node_id bigint;
    _grandchild_node_id bigint;
BEGIN
    -- Test generic tree functions
    RETURN NEXT lives_ok(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, NULL, ''
        );
    $$, 'Create root tree node');

    RETURN NEXT throws_ok(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, NULL, ''
        );
    $$, 23505, NULL, 'Create root tree node twice');

    RETURN NEXT throws_ok(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, NULL, 'test'
        );
    $$, 23514, NULL, 'Create root node with non-empty name');

    -- Get parent node id
    SELECT node_id 
    FROM projections.tree_nodes
    WHERE tree_id = 0 AND parent_id IS NULL
    INTO _parent_id;

    -- Attach child nodes to parent
    RETURN NEXT throws_ok(format(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, %s, ''
        );
    $$, _parent_id), 23514, NULL, 
    'Block creation of non-root node without name');
    
    RETURN NEXT lives_ok(format(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, %s, 'a'
        );
    $$, _parent_id), 'Attach child node to root');

    RETURN NEXT lives_ok(format(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, %s, 'b'
        );
    $$, _parent_id), 'Attach second child node to root');

    RETURN NEXT throws_ok(format(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, %s, 'a'
        );
    $$, _parent_id), 23505, NULL, 
    'Block creation of non-root node with repetitive name');

    -- Get node id for first-level node
    SELECT node_id 
    FROM projections.tree_nodes
    WHERE tree_id = 0 AND parent_id = _parent_id AND node_name = 'a'
    INTO _child_node_id;

    RETURN NEXT lives_ok(format(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, %s, 'a'
        );
    $$, _child_node_id), 'Attach second level child node');

    RETURN NEXT throws_ok(format(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, %s, 'a'
        );
    $$, _child_node_id), 23505, NULL, 
    'Block creation of second-level non-root node with repetitive name');

    RETURN NEXT results_eq(format(
    $$
        SELECT node_path 
        FROM projections.tree_nodes
        WHERE node_id = %s
    $$, _parent_id),
    $$
        SELECT '{}'::varchar[]
    $$, 'Check root node path');

    RETURN NEXT results_eq(format(
    $$
        SELECT node_path 
        FROM projections.tree_nodes
        WHERE node_id = %s
    $$, _child_node_id),
    $$
        SELECT '{""}'::varchar[]
    $$, 'Check first-level node path');

    SELECT node_id 
    FROM projections.tree_nodes
    WHERE tree_id = 0 AND parent_id = _child_node_id AND node_name = 'a'
    INTO _grandchild_node_id;

    RETURN NEXT results_eq(format(
    $$
        SELECT node_path 
        FROM projections.tree_nodes
        WHERE node_id = %s
    $$, _grandchild_node_id),
    $$
        SELECT '{"", "a"}'::varchar[]
    $$, 'Check second-level node path');

    -- Test list nodes function
    RETURN NEXT results_eq(
    $$
        SELECT COUNT(*)::integer FROM projections.list_nodes('projections.tree_nodes', 0, '{}')
    $$,
    ARRAY[1], 'Check root node path');

    RETURN NEXT results_eq(
    $$
        SELECT COUNT(*)::integer FROM projections.list_nodes('projections.tree_nodes', 0, '{""}')
    $$,
    ARRAY[2], 'Check first-level node path');

    RETURN NEXT results_eq(
    $$
        SELECT COUNT(*)::integer FROM projections.list_nodes('projections.tree_nodes', 0, '{"", "a"}')
    $$,
    ARRAY[1], 'Check first-level node path');

    RETURN NEXT results_eq(
    $$
        SELECT COUNT(*)::integer FROM projections.list_nodes('projections.tree_nodes', 0, '{"", "a", "none"}')
    $$,
    ARRAY[0], 'Check inexisting node path');

END;
$BODY$ LANGUAGE plpgsql;

-- Start testing and return the result
SELECT * FROM runtests();