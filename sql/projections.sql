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
    node_id bigserial NOT NULL,
    parent_id serial NOT NULL,
    node_name varchar,
    node_path varchar[]
);
ALTER TABLE projections.tree_nodes
    OWNER TO projections_admin;

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
    prototype_id bigserial PRIMARY KEY,
    prototype_name varchar NOT NULL,
    node_type projections.node_types NOT NULL DEFAULT 'FILE',
    uri varchar,
    -- TODO: consider context retrieval
    FOREIGN KEY (prototype_id)
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
    LIKE projections.tree_nodes,
    projection_id bigint NOT NULL,
    node_content_uri varchar NOT NULL,
    node_type projections.node_types NOT NULL DEFAULT 'FILE',
    metadata_content jsonb DEFAULT '{}',
    PRIMARY KEY (node_id),
    UNIQUE (projection_id, node_name, node_path),
    FOREIGN KEY (projection_id)
        REFERENCES projections.projections (projection_id) MATCH SIMPLE
        ON UPDATE CASCADE ON DELETE CASCADE
);
ALTER TABLE projections.projection_nodes
    OWNER TO projections_admin;
CREATE INDEX ON projections.projection_nodes(projection_id);
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
    Add node
*/
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
    _node_path varchar[];
BEGIN
    -- Create node path from parent node path
    SELECT COALESCE(array_append(node_path, node_name), ARRAY[])
    FROM projections.projection_nodes
    WHERE node_id = parent_node_id
    INTO _node_path;

    -- TODO: consider update node case
    INSERT INTO projections.projection_nodes (
        parent_id,
        node_name,
        node_path,
        projection_id,
        node_content_uri,
        node_type,
        metadata_content
    ) VALUES (
        parent_node_id,
        node_name,
        _node_path,
        projection_id,
        node_content_uri,
        node_type,
        metadata_content
    );

END;
$BODY$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION add_node() RETURNS SETOF projections.tree_nodes AS
$BODY$
DECLARE
BEGIN
END;
$BODY$ LANGUAGE plpgsql;



