/*
    Copyright 2016  Anton Bragin, Victor Svekolkin

    This file is part of Projections.

    Projections is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Projections is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Projections.  If not, see <http://www.gnu.org/licenses/>.
*/

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
    /*
        Test tree structure:
        .
        |_b
        |_a
          |_aa
    */
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
            'projections.tree_nodes', 0, %s, 'aa'
        );
    $$, _child_node_id), 'Attach second level child node');

    RETURN NEXT throws_ok(format(
    $$
        SELECT projections.add_node(
            'projections.tree_nodes', 0, %s, 'aa'
        );
    $$, _child_node_id), 23505, NULL, 
    'Block creation of second-level non-root node with repetitive name');

    RETURN NEXT results_eq(format(
    $$
        SELECT node_path IS NULL 
        FROM projections.tree_nodes
        WHERE node_id = %s
    $$, _parent_id), ARRAY[TRUE]
    , 'Check root node path');

    RETURN NEXT results_eq(format(
    $$
        SELECT node_path 
        FROM projections.tree_nodes
        WHERE node_id = %s
    $$, _child_node_id),
    $$
        SELECT '{}'::varchar[]
    $$, 'Check first-level node path');

    SELECT node_id 
    FROM projections.tree_nodes
    WHERE tree_id = 0 AND parent_id = _child_node_id AND node_name = 'aa'
    INTO _grandchild_node_id;

    RETURN NEXT results_eq(format(
    $$
        SELECT node_path 
        FROM projections.tree_nodes
        WHERE node_id = %s
    $$, _grandchild_node_id),
    $$
        SELECT '{"a"}'::varchar[]
    $$, 'Check second-level node path');

    -- Test list nodes function
    RETURN NEXT results_eq(
    $$
        SELECT COUNT(*)::integer FROM projections.list_nodes('projections.tree_nodes', 0, NULL)
    $$,
    ARRAY[1], 'Check root node path');

    RETURN NEXT results_eq(
    $$
        SELECT COUNT(*)::integer FROM projections.list_nodes('projections.tree_nodes', 0, '{}')
    $$,
    ARRAY[2], 'Check first-level node path');

    RETURN NEXT results_eq(
    $$
        SELECT COUNT(*)::integer FROM projections.list_nodes('projections.tree_nodes', 0, '{"a"}')
    $$,
    ARRAY[1], 'Check second-level node path');

    RETURN NEXT results_eq(
    $$
        SELECT COUNT(*)::integer FROM projections.list_nodes('projections.tree_nodes', 0, '{"", "a", "none"}')
    $$,
    ARRAY[0], 'Check inexisting node path');

    -- Add more nodes for testing
        /*
        Test tree structure:
        .
        |_b
        |_a
        | |_aa
        c_ca
          |_caa
          |_cab
          |_cac
          |         
          cb
    */
    SELECT projections.add_node('projections.tree_nodes', 0, _parent_id, 'c')
        INTO _child_node_id;
    SELECT projections.add_node('projections.tree_nodes', 0, _child_node_id, 'ca')
        INTO _grandchild_node_id;
    PERFORM projections.add_node('projections.tree_nodes', 0, _grandchild_node_id, 'caa');
    PERFORM projections.add_node('projections.tree_nodes', 0, _grandchild_node_id, 'cab');
    PERFORM projections.add_node('projections.tree_nodes', 0, _grandchild_node_id, 'cac');
    SELECT projections.add_node('projections.tree_nodes', 0, _child_node_id, 'cb')
        INTO _grandchild_node_id;

    RETURN NEXT results_eq(
    $$
        SELECT projections.validate_tree('projections.tree_nodes', 0)
    $$, ARRAY[TRUE], 'Check tree state');

    -- Test move node
    RETURN NEXT lives_ok(
    $$
        SELECT projections.move_node('projections.tree_nodes', 0, (
            SELECT node_id 
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = 'cab'
        ), (
            SELECT node_id
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = 'caa'
        ))
    $$, 'Check terminal node moving');

    RETURN NEXT results_eq(
    $$
        SELECT parent_id, node_path 
        FROM projections.tree_nodes
        WHERE tree_id = 0 AND node_name = 'cab'
    $$, $$
        SELECT node_id, ARRAY['c', 'ca', 'caa']::varchar[]
        FROM projections.tree_nodes
        WHERE tree_id = 0 AND node_name = 'caa'
    $$, 'Check that node and its children updated');

    RETURN NEXT lives_ok(
    $$
        SELECT projections.move_node('projections.tree_nodes', 0, (
            SELECT node_id 
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = 'ca'
        ), (
            SELECT node_id
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = ''
        ))
    $$, 'Check node moving to root');

    RETURN NEXT results_eq(
    $$
        SELECT parent_id, node_path 
        FROM projections.tree_nodes
        WHERE tree_id = 0 AND node_name = 'ca'
    $$, $$
        SELECT node_id, '{}'::varchar[]
        FROM projections.tree_nodes
        WHERE tree_id = 0 AND node_name = ''
    $$, 'Check that node and its children updated');

    RETURN NEXT results_eq(
    $$
        SELECT node_path 
        FROM projections.tree_nodes
        WHERE tree_id = 0 AND node_name = 'cab'
    $$, $$
        SELECT ARRAY['ca', 'caa']::varchar[]
    $$, 'Check that node deep children paths updated');

    RETURN NEXT lives_ok(
    $$
        SELECT projections.move_node('projections.tree_nodes', 0, (
            SELECT node_id 
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = 'ca'
        ), (
            SELECT node_id
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = 'c'
        ))
    $$, 'Check node moving from root');

    RETURN NEXT results_eq(
    $$
        SELECT node_path 
        FROM projections.tree_nodes
        WHERE tree_id = 0 AND node_name = 'cab'
    $$, $$
        SELECT ARRAY['c', 'ca', 'caa']::varchar[]
    $$, 'Check that node deep children paths updated');

    RETURN NEXT lives_ok(
    $$
        SELECT projections.move_node('projections.tree_nodes', 0, (
            SELECT node_id 
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = 'ca'
        ), (
            SELECT node_id
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = 'a'
        ))
    $$, 'Check node moving one level');

    RETURN NEXT results_eq(
    $$
        SELECT node_path 
        FROM projections.tree_nodes
        WHERE tree_id = 0 AND node_name = 'cab'
    $$, $$
        SELECT ARRAY['a', 'ca', 'caa']::varchar[]
    $$, 'Check that node deep children paths updated with one level move');

    RETURN NEXT results_eq(
    $$
        SELECT projections.validate_tree('projections.tree_nodes', 0)
    $$, ARRAY[TRUE], 'Check tree state');

    -- Test delete node
    -- Restore original structure
    PERFORM projections.move_node('projections.tree_nodes', 0, (
            SELECT node_id 
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = 'ca'
        ), (
            SELECT node_id
            FROM projections.tree_nodes
            WHERE tree_id = 0 AND node_name = 'c'
        ));
    
    RETURN NEXT lives_ok(format(
    $$
        SELECT projections.remove_node('projections.tree_nodes', 0, %s)
    $$, _child_node_id), 'Check that node can be deleted');

    RETURN NEXT results_eq(
    $$
        SELECT COUNT(*)::integer
        FROM projections.tree_nodes
        WHERE node_name LIKE 'c%'
    $$, ARRAY[0], 'Check that there are no node and its descendants');

    RETURN NEXT results_eq(
    $$
        SELECT projections.validate_tree('projections.tree_nodes', 0)
    $$, ARRAY[TRUE], 'Check tree state');

END;
$BODY$ LANGUAGE plpgsql;

-- Start testing and return the result
SELECT * FROM runtests();