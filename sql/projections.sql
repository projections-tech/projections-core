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
        -- Child nodes should be deleted by cascade
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
    projector_pid integer UNIQUE,
    prototype varchar NOT NULL -- TODO remove after prototype table added
);
ALTER TABLE projections.projections
    OWNER TO projections_admin;

COMMENT ON TABLE projections.projections IS 'Identifies Projections of resources mounted to specific path on the host system.';
COMMENT ON COLUMN projections.projections.projection_name IS 'Human-readable name used to uniquely identify projection.';
COMMENT ON COLUMN projections.projections.mount_point IS 'Absolute projection mount path.';
COMMENT ON COLUMN projections.projections.projector_pid IS 'PID of projection`s projector process.';
COMMENT ON COLUMN projections.projections.prototype IS 'Prototype of projection`s projector process.';


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
    meta_links varchar[], --NOTE this is a placehoder column, which will contain prototype id in future.
    FOREIGN KEY (tree_id)
        REFERENCES projections.projections (projection_id) MATCH SIMPLE
        ON UPDATE CASCADE ON DELETE CASCADE,
    FOREIGN KEY (tree_id, parent_id) 
        REFERENCES projections.projection_nodes (tree_id, node_id) MATCH SIMPLE
        ON DELETE CASCADE
);
ALTER TABLE projections.projection_nodes
    OWNER TO projections_admin;


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
    metadata_content jsonb,
    meta_links varchar[]
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
        metadata_content,
        meta_links
    ) = (
        add_projection_node.node_content_uri,
        add_projection_node.node_type,
        add_projection_node.metadata_content,
        add_projection_node.meta_links
    ) WHERE node_id = _projection_node_id;

    RETURN _projection_node_id;
END;
$BODY$ LANGUAGE plpgsql;

/*
This function performs path joining complaint with filesystem path, e.g. given node_name='test.bam' and
path as array {exp_1, _exp_2} it will return string '/exp_1/exp_2/test.bam'.
May be used in microcode as utility function.
*/
CREATE OR REPLACE FUNCTION join_path(path varchar[], node_name varchar)
    RETURNS varchar AS
$BODY$
DECLARE
BEGIN
    IF node_name = '' THEN
        RETURN '/';
    END IF;

    IF path = '{}' THEN
        RETURN concat('/', node_name);
    ELSE
        RETURN concat( '/', array_to_string(path, '/'), '/',node_name);
    END IF;
END;
$BODY$
LANGUAGE 'plpgsql';

COMMENT ON FUNCTION join_path(path varchar[], node_name varchar) IS
    $$Join node path with node name.

        @param: path - path of node to join.
        @param: node_name - name of a node which will be appended to path.
        @returns: joined path string.
    $$;

/*
This function returns attributes of a node on path
*/
CREATE OR REPLACE FUNCTION projections.get_projection_node_attributes(
        current_tree_id bigint,
        checked_node_path varchar
        )
    RETURNS TABLE(
        st_atime int,
        st_mtime int,
        st_ctime int,
        st_size int,
        st_mode varchar,
        st_nlink int,
        st_ino int) AS
$BODY$
BEGIN
    RETURN QUERY
    WITH node_on_path AS (
        SELECT node_id
        FROM projections.projection_nodes
        WHERE join_path(
            projections.projection_nodes.node_path,
            projections.projection_nodes.node_name ) = checked_node_path)
    SELECT projections.projection_nodes.st_atime,
        projections.projection_nodes.st_mtime,
        projections.projection_nodes.st_ctime,
        projections.projection_nodes.st_size,
        projections.projection_nodes.st_mode,
        projections.projection_nodes.st_nlink,
        projections.projection_nodes.st_ino
    FROM projections.projection_nodes, node_on_path
    WHERE projections.projection_nodes.node_id = node_on_path.node_id
    AND projections.projection_nodes.tree_id = current_tree_id;
END
$BODY$
LANGUAGE plpgsql;

COMMENT ON FUNCTION projections.get_projection_node_attributes(current_tree_id bigint, checked_node_path varchar) IS
    $$This function returns attributes of node on path.

        @param: current_tree_id - id of a projection to work with.
        @param: checked_node_path - name of a node which will be appended to path.
        @returns: list of node attributes.
    $$;

/*
This function sets attributes of a node on path
*/
CREATE OR REPLACE FUNCTION projections.set_projection_node_attributes(
    current_tree_id bigint,
    current_node_id bigint,
    current_node_size bigint,
    current_node_mode varchar
    )
    RETURNS void AS
$BODY$
DECLARE
    curr_time int := EXTRACT(epoch FROM now())::int;
BEGIN
    UPDATE projections.projection_nodes
    SET st_atime = curr_time,
        st_mtime = curr_time,
        st_ctime = curr_time,
        st_size = current_node_size,
        st_mode = current_node_mode,
        st_nlink = 0, --currently 0, subject of future changes
        st_ino = 1
    WHERE projections.projection_nodes.node_id = current_node_id
    AND projections.projection_nodes.tree_id = current_tree_id;
END;
$BODY$
LANGUAGE plpgsql;

COMMENT ON FUNCTION projections.set_projection_node_attributes(
    current_tree_id bigint,
    current_node_id bigint,
    current_node_size bigint,
    current_node_mode varchar
    ) IS
    $$This function sets attributes of node on path.

        @param: current_tree_id - id of a projection to work with.
        @param: current_node_id - id of a node which attributes will be set.
        @param: size of a node to set int.
        @param: mode of a node to set str.
        @returns: VOID.
    $$;

/*
This function performs node removal from projection_nodes table basing on node path from root,
*/
CREATE OR REPLACE FUNCTION projections.remove_node_on_path(current_tree_id bigint, node_to_remove_path varchar)
RETURNS bool AS
$BODY$
DECLARE
node_to_remove_id bigint := Null;
BEGIN
    SELECT node_id INTO node_to_remove_id
    FROM projections.projection_nodes
    WHERE join_path(projections.projection_nodes.node_path,
                    projections.projection_nodes.node_name) = node_to_remove_path AND
                    projections.projection_nodes.tree_id = current_tree_id;
    IF node_to_remove_id IS NOT NULL THEN
        PERFORM projections.remove_node('projections.projection_nodes',
                                        current_tree_id,
                                        node_to_remove_id);
        RETURN TRUE;
    ELSE
        RETURN FALSE;
    END IF;
END;
$BODY$ LANGUAGE plpgsql;

COMMENT ON FUNCTION projections.remove_node_on_path(current_tree_id bigint, node_to_remove_path varchar) IS
    $$This function removes node on given path from projection.

        @param: current_tree_id - id of a projection to work with.
        @param: node_to_remove_path - path of node to remove string.
        @returns: True if node was removed, otherwise False.
    $$;

/*
This function performs node relocation to new parent node by path
*/
CREATE OR REPLACE FUNCTION projections.move_node_on_path(
    current_tree_id bigint,
    node_to_move_path varchar,
    new_parent_path varchar
)
RETURNS bool AS
$BODY$
DECLARE
node_to_move_id bigint := Null;
new_parent_id bigint := Null;
BEGIN
    SELECT node_id INTO node_to_move_id
    FROM projections.projection_nodes
    WHERE join_path(projections.projection_nodes.node_path,
                    projections.projection_nodes.node_name) = node_to_move_path AND
                    projections.projection_nodes.tree_id = current_tree_id;

    SELECT node_id INTO new_parent_id
    FROM projections.projection_nodes
    WHERE join_path(projections.projection_nodes.node_path,
                    projections.projection_nodes.node_name) = new_parent_path AND
                    projections.projection_nodes.tree_id = current_tree_id;
    --TODO inform user which node is not present, current implementation does not allow it
    IF node_to_move_id IS NOT NULL AND new_parent_id IS NOT NULL THEN
        PERFORM projections.move_node('projections.projection_nodes',
                                      current_tree_id,
                                      node_to_move_id,
                                      new_parent_id);
        RETURN True;
    ELSE
        RETURN False;
    END IF;
END;
$BODY$ LANGUAGE plpgsql;

COMMENT ON FUNCTION projections.move_node_on_path(
    current_tree_id bigint,
    node_to_move_path varchar,
    new_parent_path varchar
) IS
    $$This function moves node on given path to new parent node on projection.

        @param: current_tree_id - id of a projection to work with.
        @param: node_to_move_path - path of node to move string.
        @param: new_parent_path - path of node new parent node string.
        @returns: bool - True if moved correctly, False otherwise.
    $$;


/*
This function reports if node on current projection tree is present in projection_nodes table
*/
CREATE OR REPLACE FUNCTION projections.projector_is_managing_path(
    current_tree_id bigint,
    current_node_path varchar
    )
  RETURNS bool AS
$BODY$
DECLARE
current_node_id bigint;
BEGIN
    IF EXISTS (
        SELECT 1
        FROM projections.projection_nodes
        WHERE join_path(node_path, node_name) = current_node_path) THEN
        RETURN True;
    ELSE
        RETURN False;
    END IF;
END;
$BODY$
LANGUAGE plpgsql;

COMMENT ON FUNCTION projections.projector_is_managing_path(
    current_tree_id bigint,
    current_node_path varchar
    ) IS
    $$This function reports if node on current projection tree is present in projection_nodes table.

        @param: current_tree_id - id of a projection to work with.
        @param: current_node_path - path of node to check string.
        @returns: bool - True if node is present in table correctly, False otherwise.
    $$;

/*
This function returns list of children nodes names on node by path
*/
CREATE OR REPLACE FUNCTION projections.get_projections_on_path(
    current_tree_id bigint,
    current_path varchar
) RETURNS SETOF varchar[] AS
$BODY$
DECLARE
node_on_path_id bigint;
BEGIN
    SELECT node_id
    INTO node_on_path_id
    FROM projections.projection_nodes
    WHERE tree_id=current_tree_id
    AND join_path(node_path, node_name) = current_path;

    RETURN QUERY SELECT array(
        SELECT node_name
        FROM projections.projection_nodes
        WHERE parent_id=node_on_path_id
        );
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.get_projections_on_path(
    current_tree_id bigint,
    current_path varchar
    ) IS
    $$This function returns list of children nodes names on node by path.

        @param: current_tree_id - id of a projection to work with.
        @param: current_path - path of parent node of projections string.
        @returns: list of projection nodes names.
    $$;


/*
This function is used to return node uri, id and st_mode for projector to use
*/
CREATE OR REPLACE FUNCTION projections.get_uri_of_projection_on_path(
    current_tree_id bigint,
    current_path varchar
) RETURNS TABLE (
    node_on_path_id bigint,
    node_on_path_uri varchar,
    node_on_path_st_mode varchar
) AS
$BODY$
DECLARE
BEGIN
    RETURN QUERY
    SELECT node_id, node_content_uri, st_mode
    FROM projections.projection_nodes
    WHERE tree_id=current_tree_id
    AND join_path(node_path, node_name) = current_path;
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.get_uri_of_projection_on_path(
    current_tree_id bigint,
    current_path varchar
    ) IS
    $$This function is used to return node uri, id and st_mode for projector to use.

        @param: current_tree_id - id of a projection to work with.
        @param: current_path - path of node to work with.
        @returns: list which contains node uri, id and st_mode.
    $$;


/*
This function lists created projections
*/
CREATE OR REPLACE FUNCTION projections.daemon_get_projections()
    RETURNS TABLE(
        _projection_id bigint,
        _projection_name varchar,
        _mount_point varchar,
        _driver varchar,
        _projector_pid int,
        _prototype varchar
) AS
$BODY$
DECLARE
BEGIN
    RETURN QUERY
    SELECT projection_id, projection_name, mount_point, driver, projector_pid, prototype FROM projections.projections;
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.daemon_get_projections() IS
    $$This function returns list of all projections in database
        @returns: _projection_id - id of projection.
        @returns: _projection_name - name of projection.
        @returns: _mount_point - projection mount point,
        @returns: _driver - name of projection`s driver.
        @returns: _projector_pid - projectors PID.
        @returns: _prototype - path to projections prototype --TODO replace with table entry
    $$;


/*
This function adds new projection into database
*/
CREATE OR REPLACE FUNCTION projections.daemon_add_projection(
    _projection_name varchar,
    _mount_point varchar,
    _driver varchar,
    _prototype varchar --temporary, until prototypes are stored in database
)    RETURNS bigint AS
$BODY$
DECLARE
_projection_id bigint;
BEGIN
    INSERT INTO projections.projections (projection_name, mount_point, driver, prototype)
    VALUES (_projection_name, _mount_point, _driver, _prototype)
    RETURNING projection_id INTO _projection_id;
    RETURN _projection_id;
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.daemon_add_projection(
    _projection_id bigint,
    _projector_pid bigint
    ) IS
    $$This function adds new projection into database.
        @param: _projection_name - name of projection.
        @param: _mount_point - projection mount point.
        @param: _driver - projection`s driver.
        @param: _prototype - projection`s prototype --TODO replace with table entry
        @returns: id of created projection.
    $$;


/*
This function sets projection`s projector PID
*/
CREATE OR REPLACE FUNCTION projections.daemon_set_projection_projector_pid(
    _projection_id bigint,
    _projector_pid bigint
)    RETURNS VOID AS
$BODY$
DECLARE
BEGIN
    UPDATE projections.projections
    SET projector_pid=_projector_pid
    WHERE projection_id = _projection_id;
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.daemon_set_projection_projector_pid(
    _projection_id bigint,
    _projector_pid bigint
    ) IS
    $$This function sets projection`s projector PID
        @param: _projection_name - name of projection which projector`s pid will be changed.
        @param: _projector_pid - PID of projector, may be NULL.
        @returns: VOID.
    $$;

/*
This function returns projection id by projection name
*/
CREATE OR REPLACE FUNCTION projections.get_projection_id_by_name(
    _projection_name varchar
)    RETURNS bigint AS
$BODY$
DECLARE
_projection_to_stop_id bigint;
BEGIN
    SELECT projection_id
    INTO _projection_to_stop_id
    FROM projections.projections
    WHERE projection_name=_projection_name;
    RETURN _projection_to_stop_id;
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.get_projection_id_by_name(
    _projection_name varchar
    ) IS
    $$This function returns projection id by projection name
        @param: _projection_name - name of projection which id will be returned.
        @returns: bigint - id of projection.
    $$;


/*
This function sets projection`s projector by given name to Null
*/
CREATE OR REPLACE FUNCTION projections.daemon_stop_projection(
    _projection_name varchar
)    RETURNS VOID AS
$BODY$
DECLARE
_projection_to_stop_id bigint;
BEGIN
    _projection_to_stop_id = projections.get_projection_id_by_name(_projection_name);
    PERFORM projections.daemon_set_projection_projector_pid(_projection_to_stop_id, NULL);
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.daemon_stop_projection(
    _projection_name varchar
    ) IS
    $$This function sets projection`s projector by given name to Null
        @param: _projection_name - name of projection to stop.
        @returns: VOID.
    $$;


/*
This function removes projection from projections table
*/
CREATE OR REPLACE FUNCTION projections.daemon_remove_projection(
    _projection_name varchar
)    RETURNS VOID AS
$BODY$
DECLARE
_projection_to_remove_id bigint;
BEGIN
    _projection_to_remove_id = projections.get_projection_id_by_name(_projection_name);
    DELETE FROM projections.projections WHERE projection_id=_projection_to_remove_id;
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.daemon_remove_projection(
    _projection_name varchar
    ) IS
    $$This function removes projection from projections table
        @param: _projection_name - name of projection to remove.
        @returns: VOID.
    $$;


/*
This function performs metadata-data binding
*/
CREATE OR REPLACE FUNCTION projections.projector_bind_metadata_to_data(
    current_tree_id bigint
) RETURNS VOID AS
$BODY$
DECLARE
_node RECORD;
_node_id bigint;
_meta_node_id bigint;
_meta_link varchar;
BEGIN
    FOR _node IN SELECT node_id, meta_links FROM projections.projection_nodes WHERE tree_id=current_tree_id LOOP
        _node_id = _node.node_id;
        -- Iterating over meta links of a node
        FOREACH _meta_link IN ARRAY _node.meta_links LOOP
            -- Executing meta_link stored in meta_links, stored query may return Null if no node were finded
            EXECUTE format(
            $$
                WITH current_node AS (
                    SELECT * FROM projections.projection_nodes WHERE node_id=$1
                )
                SELECT nodes.node_id
                FROM projections.projection_nodes AS nodes, current_node
                WHERE %s AND nodes.tree_id = $2
            $$, _meta_link)
            INTO _meta_node_id
            USING _node_id, current_tree_id;
            -- If no node finded using query continue
            IF _meta_node_id IS NOT NULL THEN
                INSERT INTO projections.projection_links (head_node_id, tail_node_id, link_name)
                VALUES (_node_id, _meta_node_id, format('link_from_%s_to_%s', _node_id, _meta_node_id));
            ELSE
                CONTINUE;
            END IF;
        END LOOP;
    END LOOP;
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.projector_bind_metadata_to_data(
    current_tree_id bigint
    ) IS
    $$This function performs metadata-data binding.
        @param: current_tree_id - id of a projection on which to perform metadata binding
        @returns: VOID.
    $$;


/*
This function performs projections search
*/
CREATE OR REPLACE FUNCTION projections.search_projections(
    current_tree_name varchar,
    search_path varchar,
    search_code varchar
) RETURNS TABLE (
    path varchar
) AS
$BODY$
DECLARE
node_on_search_path_id bigint;
current_tree_id bigint;
BEGIN
    current_tree_id = projections.get_projection_id_by_name(current_tree_name);

    SELECT node_id
    INTO node_on_search_path_id
    FROM projections.projection_nodes
    WHERE join_path(node_path, node_name) = search_path
    AND tree_id=current_tree_id;

    RETURN QUERY
    EXECUTE format(
        $$
        WITH nodes AS (
            WITH RECURSIVE descendants_ids AS (
                WITH RECURSIVE tree AS (
                    SELECT node_id, ARRAY[$1]::bigint[] AS ancestors
                    FROM projections.projection_nodes WHERE parent_id = $1

                    UNION ALL

                    SELECT projections.projection_nodes.node_id, tree.ancestors || projections.projection_nodes.parent_id
                    FROM projections.projection_nodes, tree
                    WHERE projections.projection_nodes.parent_id = tree.node_id
                )
                SELECT node_id FROM tree WHERE $1 = ANY(tree.ancestors)
            )
            SELECT projections.projection_nodes.*
            FROM projections.projection_nodes, descendants_ids
            WHERE projections.projection_nodes.node_id = descendants_ids.node_id
        ),
            links AS ( -- this table contains links for nodes on path and their attached metadata
            WITH subtree_links AS (
                SELECT head_node_id, tail_node_id, link_name
                FROM projections.projection_links, nodes
                WHERE head_node_id = nodes.node_id
                )
            SELECT
            subtree_links.head_node_id AS head_node_id,
            subtree_links.tail_node_id AS tail_node_id,
            subtree_links.link_name AS link_name,
            nodes.metadata_content AS metadata_content
            FROM subtree_links, nodes
            WHERE nodes.node_id = subtree_links.tail_node_id
        )
        SELECT join_path(nodes.node_path, nodes.node_name)
        FROM nodes
        WHERE %s;
        $$, search_code)
    USING node_on_search_path_id;
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.search_projections(
    current_tree_name varchar,
    search_path varchar,
    search_code varchar
    ) IS
    $$This function is used to return node uri, id and st_mode for projector to use.
        @param: current_tree_name - name of a projections tree on which to performs search.
        @param: search_path - path of parent node which children nodes will be filtered.
        @search_code: search query SQL code.
        @returns: VOID.
    $$;


/*
This function checks if mount point is already occupied by another projection
*/
CREATE OR REPLACE FUNCTION projections.check_mount_point(
    _mount_point varchar
)    RETURNS bool AS
$BODY$
DECLARE
_is_mount_point_occupied bool;
BEGIN
    SELECT EXISTS (
        SELECT 1
        FROM projections.projections
        WHERE mount_point = _mount_point
        )
    INTO _is_mount_point_occupied;
    RETURN _is_mount_point_occupied;
END;
$BODY$ LANGUAGE 'plpgsql';

COMMENT ON FUNCTION projections.check_mount_point(
    _mount_point varchar
    ) IS
    $$This function checks if mount point is already occupied by another projection
        @param: _mount_point - mount point path.
        @returns: boll - True if occupied, False otherwise.
    $$;