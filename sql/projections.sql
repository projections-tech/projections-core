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

/*
    This script initiates database required Projections system to run.
    It assumes PostgreSQL 9.4+ and should be run on clean database.
    
    IMPORTANT: this scripts will delete all Projections data from database if any present!
*/

-- Drop database schema
DROP SCHEMA IF EXISTS projections CASCADE;

-- Drop and create main database role
DROP ROLE IF EXISTS projections_admin;
CREATE ROLE projections_admin;

CREATE SCHEMA projections 
    AUTHORIZATION projections_admin;

/*
    Tables and functions for generic tree operations.
*/
-- Generic tree tables
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
    CHECK (node_id != parent_id),
    FOREIGN KEY (tree_id, parent_id) 
        REFERENCES projections.tree_nodes (tree_id, node_id) MATCH SIMPLE
        ON DELETE CASCADE
);
ALTER TABLE projections.tree_nodes
    OWNER TO projections_admin;
-- This index prevents from having multiply root nodes for one tree
CREATE UNIQUE INDEX ON projections.tree_nodes (tree_id, node_name) 
    WHERE node_path IS NULL;

CREATE INDEX ON projections.tree_nodes (tree_id);

COMMENT ON TABLE projections.tree_nodes IS 
    $$Generic table describing tree nodes. 

        Designed as 'abstract' node class for extending by prototype and projection nodes.
        Records of this type are used by generic tree functions.

        LIKE projections.tree_nodes INCLUDING CONSTRAINTS INCLUDING INDEXES INCLUDING DEFAULTS are expected to be used
        when inheriting from the table.
    $$;

-- TODO: add links and linking / search operations
    
-- Generic tree functions
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

COMMENT ON FUNCTION projections.add_node(table_name varchar, tree_id bigint, parent_id bigint, node_name varchar) IS 
    $$Adds new node to the specified tree with the parent provided.

        @param: table_name - represents string name of the table LIKE projections.tree_nodes to perform operation on.
        @param: tree_id - id of the tree to insert node to.
        @param: parent_id - node_id of parent node. Should be NULL for root node.
        @param: node_name - name of node to create. This name will be included in paths of all child nodes. Should be '' for root
            and not empty for any non-root node.
        @returns: node_id of created node.
    $$;


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

COMMENT ON FUNCTION projections.list_nodes(table_name varchar, tree_id bigint, node_path varchar[]) IS 
    $$List nodes  on the path specified.

        @param: table_name - represents string name of the table LIKE projections.tree_nodes to perform operation on.
        @param: tree_id - id of the tree to insert node to.
        @param: node_path – path to search nodes on.
        @returns: list of projections.tree_nodes objects that may be joined with table_name records.
    $$;


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


CREATE OR REPLACE FUNCTION projections.aggregate_node_names(
    table_name varchar,
    tree_id bigint
) RETURNS TABLE (
    node_id bigint,
    node_name varchar,
    node_path varchar[],
    node_names varchar[],
    level integer
) AS
$BODY$
DECLARE
BEGIN
    RETURN QUERY EXECUTE format($$
        (WITH RECURSIVE parent AS (
            SELECT node_id, node_name, node_path, ARRAY[]::varchar[] AS node_names, 0 AS level
            FROM %1$s
            WHERE tree_id = $1 AND parent_id IS NULL
        UNION ALL
            SELECT 
                target.node_id,
                target.node_name,
                target.node_path,
                parent.node_names || target.node_name AS node_names,
                parent.level + 1 AS level
            FROM %1$s AS target
                INNER JOIN parent ON (target.parent_id = parent.node_id)
        ) SELECT * FROM parent)
    $$, table_name)
    USING tree_id;
END;        
$BODY$ LANGUAGE plpgsql;

COMMENT ON FUNCTION projections.aggregate_node_names(table_name varchar, tree_id bigint) IS 
    $$Aggregate node names from root node to all descendants.

        May be used to test whether the tree specified by id is in consistent state.

        @param: table_name - represents string name of the table LIKE projections.tree_nodes to perform operation on.
        @param: tree_id - id of the tree to insert node to.
        @returns: node_path - path of node
        @returns: node_names - aggragate of node names from root to current node
        @returns: level - number of node parents on the way to root
    $$;


CREATE OR REPLACE FUNCTION projections.validate_tree(
    table_name varchar,
    tree_id bigint,
    OUT is_valid boolean
) AS
$BODY$
DECLARE
BEGIN
    -- Check path for all nodes except the root 
    WITH comparison AS (
        SELECT node_path || node_name = node_names AS path_corresp
        FROM projections.aggregate_node_names(table_name, tree_id)
        WHERE node_path IS NOT NULL
    ) SELECT every(comparison.path_corresp)
    FROM comparison
    INTO is_valid;
END;
$BODY$ LANGUAGE plpgsql;

COMMENT ON FUNCTION projections.validate_tree(table_name varchar, tree_id bigint) IS 
    $$Check whether node paths comsistent with tree hierarchy for the tree provided.
    
        Tests whether the tree specified by id is in consistent state.

        @param: table_name - represents string name of the table LIKE projections.tree_nodes to perform operation on.
        @param: tree_id - id of the tree to insert node to.
        @returns: boolean - true if tree is valid, false otherwise.
    $$;

/*
    Tables and functions for projection-specific tables and functions definitions.
*/

/*
    Filesystem node types. See: 
        http://man7.org/linux/man-pages/man2/stat.2.html 
    for additional details.
    Only 'REG' and 'DIR' are supported at the moment.
*/ 
CREATE TYPE projections.node_types AS ENUM (
    'SOCK',
    'LNK',
    'REG',
    'BLK',
    'DIR',
    'CHR',
    'FIFO'
);

-- Tables for prototypes data storage
-- TODO: do prototype tables revision
CREATE TABLE projections.prototypes (
    prototype_id bigserial PRIMARY KEY,
    prototype_name varchar NOT NULL
);
ALTER TABLE projections.prototypes
    OWNER TO projections_admin;

CREATE TABLE projections.prototype_nodes (
    LIKE projections.tree_nodes,
    prototype_name varchar NOT NULL,
    node_type projections.node_types NOT NULL DEFAULT 'REG',
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
    LIKE projections.tree_nodes INCLUDING CONSTRAINTS INCLUDING INDEXES INCLUDING DEFAULTS,
    node_content_uri varchar,
    node_type projections.node_types DEFAULT 'REG',
    st_atime int,
    st_mtime int,
    st_ctime int,
    st_size int DEFAULT 1,
    st_mode varchar,
    st_nlink int,
    st_ino int,
    metadata_content jsonb DEFAULT '{}',
    FOREIGN KEY (tree_id)
        REFERENCES projections.projections (projection_id) MATCH SIMPLE
        ON UPDATE CASCADE ON DELETE CASCADE,
    FOREIGN KEY (tree_id, parent_id) 
        REFERENCES projections.projection_nodes (tree_id, node_id) MATCH SIMPLE
        ON DELETE CASCADE
);
ALTER TABLE projections.projection_nodes
    OWNER TO projections_admin;

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

-- Functions for projection tree operation
/*
    This function is just a synthetic sugar.
    projections.add_node may be used for tree operations while fields setup may be performed via updates.
*/
CREATE OR REPLACE FUNCTION projections.add_projection_node(
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
        add_projection_node.node_content_uri,
        add_projection_node.node_type,
        add_projection_node.metadata_content
    ) WHERE node_id = _projection_node_id;

    RETURN _projection_node_id;
END;
$BODY$ LANGUAGE plpgsql;
