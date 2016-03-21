-- Tables creation

CREATE TABLE projections_table (
                projection_name varchar PRIMARY KEY UNIQUE NOT NULL,
                mount_path varchar UNIQUE,
                projector_pid int UNIQUE
                );

COMMENT ON TABLE projections_table IS 'This table stores projections names mount points and PID`s of their projectors.';
COMMENT ON COLUMN projections_table.projection_name IS 'Name of projection, must be unique.';
COMMENT ON COLUMN projections_table.mount_path IS 'Absolute path on which projection is mounted.';
COMMENT ON COLUMN projections_table.projector_pid IS 'PID of projection`s projector.';


CREATE TABLE tree_table (
                node_id serial PRIMARY KEY,
                projection_name varchar REFERENCES projections_table(projection_name) ON DELETE CASCADE NOT NULL,
                parent_id integer,
                name varchar NOT NULL,
                uri varchar NOT NULL,
                path varchar[][] NOT NULL,
                type varchar NOT NULL,
                meta_links varchar[][]
                );

COMMENT ON TABLE tree_table IS 'This table holds projection tree structure.';
COMMENT ON COLUMN tree_table.node_id IS 'Id integer of projection node.';
COMMENT ON COLUMN tree_table.projection_name IS 'Name string of projection.';
COMMENT ON COLUMN tree_table.parent_id IS 'Id integer of node`s parent node.';
COMMENT ON COLUMN tree_table.name IS 'Name string of node.';
COMMENT ON COLUMN tree_table.uri IS 'Node resource uri string.';
COMMENT ON COLUMN tree_table.path IS 'Path to node list of strings.';
COMMENT ON COLUMN tree_table.type IS 'Node type string dir or file';
COMMENT ON COLUMN tree_table.meta_links IS 'List of metadata links strings.';

CREATE TABLE metadata_table (
                meta_id serial PRIMARY KEY,
                node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE NOT NULL,
                parent_node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE NOT NULL,
                meta_contents jsonb NOT NULL
                );

COMMENT ON TABLE metadata_table IS 'This table holds data-metadata binding relation and metadata contents.';
COMMENT ON COLUMN metadata_table.meta_id IS 'Id integer of data-metadata relation.';
COMMENT ON COLUMN metadata_table.node_id IS 'Id integer of metadata node.';
COMMENT ON COLUMN metadata_table.parent_node_id IS 'Id of data node to which metadata is binded.';
COMMENT ON COLUMN metadata_table.meta_contents IS 'Metadata node contents as jsonb.';

-- TODO change name to attributes_table after enough functions converted to SQL
CREATE TABLE projection_attributes_table (
                attr_id serial PRIMARY KEY,
                node_id integer REFERENCES tree_table(node_id) ON DELETE CASCADE,
                st_atime int NOT NULL,
                st_mtime int NOT NULL,
                st_ctime int NOT NULL,
                st_size int DEFAULT 1 NOT NULL,
                st_mode varchar NOT NULL,
                st_nlink int NOT NULL,
                st_ino int Not NULL
                );

COMMENT ON TABLE projection_attributes_table IS 'This table holds projection nodes filesystem attributes.';
COMMENT ON COLUMN projection_attributes_table.attr_id IS 'Id int of attribute.';
COMMENT ON COLUMN projection_attributes_table.node_id IS 'Node which attributes are stored id int.';
COMMENT ON COLUMN projection_attributes_table.st_atime IS 'Time of last access.';
COMMENT ON COLUMN projection_attributes_table.st_mtime IS 'Time of last modification.';
COMMENT ON COLUMN projection_attributes_table.st_ctime IS 'Time of last status change.';
COMMENT ON COLUMN projection_attributes_table.st_size IS 'Size of node.';
COMMENT ON COLUMN projection_attributes_table.st_nlink IS 'Number of hard links.';
COMMENT ON COLUMN projection_attributes_table.st_ino IS 'Inode number.';

-- Defining functions
