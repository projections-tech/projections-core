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
    -- Add test projection
    INSERT INTO projections.projections (
        projection_name, mount_point, driver
    ) VALUES (
        'test_projection', '/home/test', 'fs_driver'
    );
	
END;
$$ LANGUAGE plpgsql;

-- -- Main test function
CREATE OR REPLACE FUNCTION test_projections_functions() RETURNS SETOF text AS 
$BODY$
DECLARE
    _projection_id bigint;
    _root_node_id bigint;
    _child_node_id bigint;
    _second_child_node_id bigint;
BEGIN
    -- Get projection id for testing
    SELECT projection_id 
    FROM projections.projections
    WHERE projection_name = 'test_projection'
    INTO _projection_id;

    RETURN NEXT throws_ok(format($$
        RETURN NEXT (SELECT projections.add_projection_node(_projection_id, 'a', 1, 'http://projections.tech', 'DIR', NULL));
    $$, _projection_id), NULL, NULL, 'Create projection node before root node creation');

    RETURN NEXT lives_ok(format($$
        SELECT projections.add_projection_node(%s, '', NULL, 'http://projections.tech/', 'DIR', NULL, '{}');
    $$, _projection_id), 'Create projection root node');

    SELECT node_id 
    FROM projections.projection_nodes
    WHERE tree_id = _projection_id AND parent_id IS NULL
    INTO _root_node_id;

    RETURN NEXT lives_ok(format($$
        SELECT projections.add_projection_node(%s, 'a_dir', %s, 'http://projections.tech/a_dir', 'DIR', NULL, '{}');
    $$, _projection_id, _root_node_id), 'Create projection first-level node');

    SELECT node_id 
    FROM projections.projection_nodes
    WHERE tree_id = _projection_id AND node_name = 'a_dir'
    INTO _child_node_id;

    RETURN NEXT lives_ok(format($$
        SELECT projections.add_projection_node(%s, 'aa_dir', %s, 'http://projections.tech/a_dir/aa_dir', 'DIR', NULL, '{}');
    $$, _projection_id, _child_node_id), 'Create projection second-level node');

    RETURN NEXT lives_ok(format($$
        SELECT projections.add_node('projections.projection_nodes', %s, %s, 'ab_dir');
    $$, _projection_id, _child_node_id), 'Create projection second-level node by generic function');

    -- Add more nodes
    SELECT projections.add_projection_node(_projection_id, 'b_dir', _root_node_id, 'http://projections.tech/b_dir', 'DIR', NULL, '{}')
        INTO _second_child_node_id;
    PERFORM projections.add_projection_node(_projection_id, 'c', _root_node_id, 'http://projections.tech/c', 'REG', NULL, '{}');
    PERFORM projections.add_projection_node(_projection_id, 'aa', _child_node_id, 'http://projections.tech/a_dir/aa', 'REG', NULL, '{}');

    RETURN NEXT results_eq(format($$
        SELECT projections.validate_tree('projections.projection_nodes', %s)
    $$, _projection_id), ARRAY[TRUE], 'Check tree state');

    RETURN NEXT lives_ok(format($$
        SELECT projections.move_node('projections.projection_nodes', %s, %s, %s);
    $$, _projection_id, _child_node_id, _second_child_node_id), 'Move projection tree node');

    RETURN NEXT results_eq(format($$
        SELECT projections.validate_tree('projections.projection_nodes', %s)
    $$, _projection_id), ARRAY[TRUE], 'Check tree state');

    RETURN NEXT lives_ok(format($$
        SELECT projections.move_node('projections.projection_nodes', %s, %s, %s);
    $$, _projection_id, _child_node_id, _second_child_node_id), 'Move projection tree node');

    RETURN NEXT results_eq(format($$
        SELECT projections.validate_tree('projections.projection_nodes', %s)
    $$, _projection_id), ARRAY[TRUE], 'Check tree state');



END;
$BODY$ LANGUAGE plpgsql;

-- Start testing and return the result
SELECT * FROM runtests('projections_functions');