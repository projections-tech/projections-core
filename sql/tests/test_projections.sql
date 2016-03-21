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
    _projection_id bigint;
BEGIN
    RETURN NEXT lives_ok(
    $$
        INSERT INTO projections.projections (
            projection_name,
            mount_point,
            driver
        ) VALUES (
            'test_projection',
            '/home/test/',
            'fs_driver'
        );
    $$, 'Create projection');

    -- Get newly created projection ID
    _projection_id := (
        SELECT projection_id 
        FROM projections.projections
        WHERE projection_name = 'test_projection'
    );

	RETURN;
END;
$BODY$ LANGUAGE plpgsql;

-- Start testing and return the result
SELECT * FROM runtests();