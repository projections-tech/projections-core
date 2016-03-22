-- IMPORTANT!: this scripts will delete all data from database if any present!

-- Drop database schema
DROP SCHEMA IF EXISTS projections CASCADE;

-- Drop and create main database role
DROP ROLE IF EXISTS projections_admin;
CREATE ROLE projections_admin;

CREATE SCHEMA projections 
    AUTHORIZATION projections_admin;

-- Generic nodes table
CREATE TABLE projections.tree_nodes (
    node_id bigserial PRIMARY KEY,
    tree_id bigint NOT NULL,
    parent_id bigint,
    node_name varchar NOT NULL,
    node_path varchar[],
    UNIQUE (tree_id, node_path, node_name),
    UNIQUE (tree_id, node_id),
    CHECK (node_name != '' OR parent_id IS NULL),
    CHECK (node_name = '' OR parent_id IS NOT NULL),
    FOREIGN KEY (tree_id, parent_id) 
        REFERENCES projections.tree_nodes (tree_id, node_id) MATCH SIMPLE
        ON DELETE CASCADE
);
ALTER TABLE projections.tree_nodes
    OWNER TO projections_admin;
-- This index prevents from having multiply root nodes for one tree
CREATE UNIQUE INDEX ON projections.tree_nodes (tree_id, node_name) 
    WHERE node_path IS NULL; 

-- TODO: harmonize with POSIX filesystem object types
CREATE TYPE projections.node_types AS ENUM (
    'FILE',
    'FOLDER'
);

-- Tables for prototypes data storage
CREATE TABLE projections.prototypes (
    prototype_id bigserial PRIMARY KEY,
    prototype_name varchar NOT NULL
);
ALTER TABLE projections.prototypes
    OWNER TO projections_admin;

CREATE TABLE projections.prototype_nodes (
    LIKE projections.tree_nodes,
    prototype_name varchar NOT NULL,
    node_type projections.node_types NOT NULL DEFAULT 'FILE',
    uri varchar,
    -- TODO: consider context retrieval
    PRIMARY KEY (node_id),
    FOREIGN KEY (tree_id)
        REFERENCES projections.prototypes (prototype_id) MATCH SIMPLE
        ON UPDATE CASCADE ON DELETE CASCADE
);
ALTER TABLE projections.prototype_nodes
    OWNER TO projections_admin;  

-- Projection tables
CREATE TABLE projections.projections (
    projection_id bigserial PRIMARY KEY,
    projection_name varchar NOT NULL UNIQUE,
    mount_point varchar NOT NULL UNIQUE,
    driver varchar NOT NULL,    -- May reference additional drivers table records
    projector_pid integer UNIQUE
);
ALTER TABLE projections.projections
    OWNER TO projections_admin;

COMMENT ON TABLE projections_table IS 'Identifies Projections of resources mounted to specific path on the host system.';
COMMENT ON COLUMN projections_table.projection_name IS 'Human-readable name used to uniquely identify projection.';
COMMENT ON COLUMN projections_table.mount_path IS 'Absolute projection mount path.';
COMMENT ON COLUMN projections_table.projector_pid IS 'PID of projection`s projector process.';

CREATE TABLE projections.projection_nodes (
    LIKE projections.tree_nodes INCLUDING CONSTRAINTS INCLUDING INDEXES,
    node_content_uri varchar,
    node_type projections.node_types DEFAULT 'FILE',
    metadata_content jsonb DEFAULT '{}',
    FOREIGN KEY (tree_id)
        REFERENCES projections.projections (projection_id) MATCH SIMPLE
        ON UPDATE CASCADE ON DELETE CASCADE
);
ALTER TABLE projections.projection_nodes
    OWNER TO projections_admin;
CREATE INDEX ON projections.projection_nodes(tree_id);
-- TODO: Add indexes required for tree operations
COMMENT ON TABLE tree_table IS 'Holds projection node records each reflecting some data object projected as file or folder.';

CREATE TABLE projections.projection_links (
    projection_link_id bigserial PRIMARY KEY,
    head_node_id bigint NOT NULL,
    tail_node_id bigint NOT NULL,
    link_name varchar NOT NULL,
    UNIQUE (head_node_id, tail_node_id, link_name),
    FOREIGN KEY (head_node_id)
        REFERENCES projections.projection_nodes (node_id) MATCH SIMPLE
        ON UPDATE CASCADE ON DELETE CASCADE,
    FOREIGN KEY (tail_node_id)
        REFERENCES projections.projection_nodes (node_id) MATCH SIMPLE
        ON UPDATE CASCADE ON DELETE CASCADE
    -- NOTE: with this design it is in principle possible to link objects from different projections.
    -- Review this and restrict if needed.
);
ALTER TABLE projections.projection_links
    OWNER TO projections_admin;
COMMENT ON TABLE projections.projection_links IS 'Holds named and directed links between different projection objects (such as metadata links).';

-- Consider merging this data with projection_nodes since the tables are connected with one-to-one relationship 
CREATE TABLE projections.projection_node_fs_attributes (
    projection_node_fs_properties_id bigserial PRIMARY KEY,
    node_id integer NOT NULL UNIQUE,
    st_atime int NOT NULL,
    st_mtime int NOT NULL,
    st_ctime int NOT NULL,
    st_size int DEFAULT 1 NOT NULL,
    st_mode varchar NOT NULL,
    st_nlink int NOT NULL,
    st_ino int NOT NULL,
    FOREIGN KEY (node_id)
        REFERENCES projections.projection_nodes (node_id) MATCH SIMPLE
        ON UPDATE CASCADE ON DELETE CASCADE
);
ALTER TABLE projections.projection_node_fs_attributes
    OWNER TO projections_admin;

COMMENT ON TABLE projections.projection_node_fs_attributes IS 'Holds projection node filesystem attributes.';

-- Generic tree functions
/*
    Add node to the table specified and return node_id.
    It is expected that table LIKE projections.tree_nodes
*/
CREATE OR REPLACE FUNCTION projections.add_node(
    table_name varchar,
    tree_id bigint,
    parent_id bigint,
    node_name varchar
) RETURNS bigint AS
$BODY$
DECLARE
    _node_path varchar[];
    _node_id bigint;
BEGIN
    -- Create node path from parent node path
    EXECUTE format($$
        SELECT 
            CASE WHEN node_name IS NULL OR node_name = ''
                THEN
                    '{}'
                ELSE
                    array_append(node_path, node_name)
            END
        FROM %s
        WHERE node_id = $1
    $$, table_name)
    INTO _node_path
    USING parent_id;

    EXECUTE format($$
        INSERT INTO %s (
            tree_id,
            parent_id,
            node_name,
            node_path
        ) VALUES (
            $1,
            $2,
            $3,
            $4
        ) RETURNING node_id
    $$, table_name)
    INTO _node_id
    USING tree_id, parent_id, node_name, _node_path;
    -- Return id of node created
    RETURN _node_id;
END;
$BODY$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION projections.list_nodes(
    table_name varchar,
    tree_id bigint,
    node_path varchar[]
) RETURNS SETOF projections.tree_nodes AS
$BODY$
DECLARE
BEGIN
    CASE WHEN node_path IS NOT NULL THEN
        RETURN QUERY
            EXECUTE format($$
                SELECT node_id, tree_id, parent_id, node_name, node_path
                FROM %s
                WHERE tree_id = $1 AND node_path = $2
            $$, table_name)
            USING tree_id, node_path;
    ELSE
        RETURN QUERY
            EXECUTE format($$
                SELECT node_id, tree_id, parent_id, node_name, node_path
                FROM %s
                WHERE tree_id = $1 AND node_path IS NULL
            $$, table_name)
            USING tree_id;
    END CASE;
END;
$BODY$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION projections.move_node(
    table_name varchar,
    tree_id bigint,
    node_id bigint,
    new_parent_id bigint
) RETURNS SETOF projections.tree_nodes AS
$BODY$
DECLARE
    _path varchar[];
    _path_length integer;
    _new_path varchar[];
BEGIN
    -- Get path and path length of current node
    EXECUTE format($$
        SELECT node_path || node_name, cardinality(node_path) + 1 
        FROM %s
        WHERE node_id = $1
    $$, table_name)
    INTO _path, _path_length
    USING node_id;

    -- Get path and length of new parent node
    EXECUTE format($$
        SELECT 
            CASE node_name 
            WHEN '' THEN
                -- Moving to root
                '{}'
            ELSE
                node_path || node_name
            END
        FROM %s
        WHERE node_id = $1
    $$, table_name)
    INTO _new_path
    USING new_parent_id;

    RAISE NOTICE 'new id: %; new path: %', new_parent_id, _new_path;   

    -- Change parent id and path of current node
    EXECUTE format($$
            UPDATE %s SET (
                parent_id,
                node_path
            ) = (
                $1,
                $2
            ) WHERE node_id = $3
    $$, table_name)
    USING new_parent_id, _new_path, node_id;

    -- Update descendants paths
    EXECUTE format($$
            UPDATE %s SET node_path = $1 || node_path[$3:cardinality(node_path)]
            WHERE tree_id = $2 AND node_path[1:$3] = $4
    $$, table_name)
    USING _new_path, tree_id, _path_length, _path;

    RETURN QUERY
        EXECUTE format($$
            SELECT node_id, tree_id, parent_id, node_name, node_path
            FROM %s
            WHERE tree_id = $1 AND node_id = $2
        $$, table_name)
        USING tree_id, node_id;
END;
$BODY$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION projections.rename_node(
    table_name varchar,
    tree_id bigint,
    node_id bigint,
    new_name varchar
) RETURNS VOID AS
$BODY$
DECLARE
BEGIN
    -- TODO: add implementation!
END;
$BODY$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION projections.remove_node(
    table_name varchar,
    tree_id bigint,
    node_id bigint
) RETURNS VOID AS
$BODY$
DECLARE
BEGIN
        EXECUTE format($$
            DELETE FROM %s
            WHERE tree_id = $1 AND node_id = $2
        $$, table_name)
        USING tree_id, node_id;
        -- Child nodes should be deleted by cascase
END;
$BODY$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION projections.validate_tree(
    table_name varchar,
    tree_id bigint
) RETURNS boolean AS
$BODY$
DECLARE
BEGIN
    /*
        TODO: add implementation!
    */
END;
$BODY$ LANGUAGE plpgsql;


-- Functions for projection tree operation
CREATE OR REPLACE FUNCTION add_projection_node(
    projection_id bigint,
    node_name varchar,
    parent_node_id bigint,
    node_content_uri varchar,
    node_type projections.node_types,
    metadata_content jsonb
    ) RETURNS bigint AS
$BODY$
DECLARE
    _projection_node_id bigint;
BEGIN
    -- Create node
    SELECT projections.add_node(
        'projections.projection_nodes',
        projection_id,
        parent_node_id,
        node_name)
    INTO _projection_node_id;

    UPDATE projections.projection_nodes SET (
        node_content_uri,
        node_type,
        metadata_content
    ) = (
        node_content_uri,
        node_type,
        metadata_content
    ) WHERE node_id = _projection_node_id;

END;
$BODY$ LANGUAGE plpgsql;

