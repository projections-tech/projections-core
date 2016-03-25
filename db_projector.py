#    Copyright 2016  Anton Bragin, Victor Svekolkin
#
#    This file is part of Projections.
#
#    Projections is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Projections is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Projections.  If not, see <http://www.gnu.org/licenses/>.

import copy
import getpass
import io
import logging
import logging.config
import os
import stat
import types

import objectpath
import psycopg2
import psycopg2.extras

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('db_projector')


class DBProjector:
    """
    This class contains test interface to PostgresSQL database, which will be used to store and work with projections.
    Current implementation is subject of future rewrites and optimisations,
    """

    def __init__(self, projection_id, projection_driver, prototypes_tree, root_uri):
        """
        This method initializes database which will hold projection tree and associated metadata
        :param projection_id: name of projection string
        :param projection_driver: projection driver instance
        :param prototypes_tree: prototypes tree instance
        :param root_uri: root uri of a projection
        """
        self.projection_id = projection_id
        self.prototypes_tree = prototypes_tree
        self.root_uri = root_uri

        # Initializing projection driver
        self.projection_driver = projection_driver

        # Opening connection with database
        self.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))

        # Creating cursor, which will be used to interact with database
        self.cursor = self.db_connection.cursor()

        self.build_projection()
        logger.debug('Initialized projection: {}'.format(self.projection_id))

    def __del__(self):
        """
        This method closes cursor and database connection at call to destructor method
        """
        self.cursor.close()
        self.db_connection.close()

    def __split_path(self, path):
        """
        Split path string in parts preserving path root
        :param path: path string
        :return: list of path parts strings
        """
        temp = []
        while True:
            head, tail = os.path.split(path)
            if head == '/' or head == '':
                if tail == '':
                    return [head]
                else:
                    temp.append(tail)
                    if not head == '':
                        temp.append(head)
                    return temp[::-1]
            temp.append(tail)
            path = head

    def build_projection(self):
        """
        This method builds projection
        """
        # Add root node to tables
        self.cursor.execute("""

        SELECT projections.add_projection_node(%(projection_id)s, %(node_name)s, %(parent_id)s, %(uri)s, %(type)s, %(metadata_content)s)

        """, {'projection_id': self.projection_id,
              'parent_id': None,
              'node_name': '',
              'type': 'DIR',
              'uri': self.root_uri,
              'metadata_content': None})

        # After insertion cursor returns id of root node, which will be used in tree building
        new_parent_id = self.cursor.fetchone()

        self.db_connection.commit()

        # Starting projections tree construction
        self.db_build_tree({'/': self.prototypes_tree}, root_projection_uri=self.root_uri, parent_id=new_parent_id)
        # self.bind_metadata_to_data()

    def db_build_tree(self, prototypes, parent_id=None, root_projection_uri=None):
        """
        This method recursively builds projections file structure in tree_table from prototype tree
        :param prototypes: Prototype Tree object
        :param parent_id: node_id of parent projection in tree_table
        :param projection_path: path to current projection string
        :param root_projection_uri: root projection uri string
        :param parent_id_for_meta: id of parent node for projections with metadata
        """
        environment = None
        content = None
        path = os.path

        # This is environment in which projections are created (parent_projection content)
        # TODO: in many cases it means double request to parent projection resource so it should be optimized
        # We don`t want to change driver contents, hence we made deep copy of dict
        environment = copy.deepcopy(self.projection_driver.get_uri_contents_as_dict(root_projection_uri))

        # For every prototype in collection try to create corresponding projections
        for key, prototype in prototypes.items():
            # Set current prototype context to current environment for children node to use
            prototype.context = environment
            # Get context of current node from contexts of parent nodes
            context = prototype.get_context()
            context = context[::-1]

            # Adding context of upper level prototypes for lower level projections to use
            environment['context'] = context

            # Creating tree of environment contents which will be parsed by ObjectPath
            tree = objectpath.Tree(environment)

            URIs = tree.execute(prototype.uri)

            # Object path sometimes returns generator if user uses selectors, for consistency expand it using
            # list comprehension
            if isinstance(URIs, types.GeneratorType):
                URIs = [el for el in URIs]
            # Treating URIs as list for consistency
            if not isinstance(URIs, list):
                URIs = [URIs]
            logger.debug('Projection uri: %s', URIs)
            # We get projection URIs based on environment and prototype properties
            # Every URI corresponds to projection object
            for uri in URIs:
                # Get content for a projection
                # We don`t want to change driver contents, hence we made deep copy
                content = copy.deepcopy(self.projection_driver.get_uri_contents_as_dict(uri))

                # Adding environment to use by prototype
                content['environment'] = environment
                content['context'] = context

                # Creating tree which will be parsed by ObjectPath
                tree = objectpath.Tree(content)
                name = tree.execute(prototype.name)
                logger.debug('Projection name: %s', name)

                # Object path sometimes returns generator if user uses selectors, for consistency expand it using
                # list comprehension
                if isinstance(name, types.GeneratorType):
                    name = [el for el in name]

                if prototype.type == 'transparent':
                    # If prototype is transparent skip projection building and pass parent node id as current level id
                    self.db_build_tree(prototype.children, parent_id, root_projection_uri=uri)
                    continue
                else:
                    current_parent_id = parent_id

                prototype_types_binding = {
                    'file': 'REG',
                    'directory': 'DIR'
                }

                # Inserting node into projections tree on completion this command returns inserted node id
                # which will be passed to lower level nodes
                self.cursor.execute("""

                SELECT projections.add_projection_node(%(projection_id)s, %(node_name)s, %(parent_id)s, %(uri)s, %(type)s, %(metadata_content)s)

                """, {'projection_id': self.projection_id,
                      'parent_id': current_parent_id,
                      'node_name': name,
                      'type': prototype_types_binding[prototype.type],
                      'uri': self.root_uri,
                      'metadata_content': None})

                # Fetching inserted projection id which will be parent id for lower level projections
                new_parent_id = self.cursor.fetchone()

                # If new_parent_id is None, when projection already exists in table
                if new_parent_id is not None:
                    self.db_build_tree(prototype.children, new_parent_id, root_projection_uri=uri)
                    # Commit all changes
                    self.db_connection.commit()

    def bind_metadata_to_data(self):
        """
        Performs data-metadata binding according to meta_link
        """
        # List all projections, getting their id`s and lists of meta links
        self.cursor.execute("""
        SELECT node_id, meta_links FROM tree_table WHERE projection_name = %s
        """, (self.projection_id,))

        projections_list = [row for row in self.cursor]

        # Perform metadata binding for each projection according to commands in meta_links
        for node_id, meta_links in projections_list:
            for meta_link in meta_links:
                if meta_link is not None:

                    self.cursor.execute("""
                    WITH meta_id AS (
                        WITH current_node AS (
                            SELECT * FROM tree_table WHERE node_id=%(current_node_id)s
                        )

                        SELECT tree_table.node_id
                        FROM tree_table, current_node
                        WHERE {meta_link} AND tree_table.projection_name = %(projection_name)s
                    )
                    INSERT INTO metadata_table (parent_node_id, node_id) SELECT %(current_node_id)s, meta_id.node_id
                    FROM meta_id
                    WHERE NOT EXISTS (SELECT 1 FROM metadata_table, meta_id
                    WHERE metadata_table.node_id = meta_id.node_id AND metadata_table.parent_node_id = %(current_node_id)s)
                    RETURNING meta_id
                    """.format(meta_link=meta_link, ), {'current_node_id': node_id,
                                                        'projection_name': self.projection_id})

                    meta_node_id = self.cursor.fetchone()
                    # Adding meta contents jsonb to projection if insertion were performed
                    if meta_node_id is not None:
                        self.cursor.execute("""
                        SELECT uri FROM tree_table WHERE node_id=(SELECT node_id FROM metadata_table WHERE meta_id = %s)
                        """, (meta_node_id,))

                        meta_node_uri = self.cursor.fetchone()[0]

                        metadata_contents = self.projection_driver.get_uri_contents_as_dict(meta_node_uri)

                        self.cursor.execute("""
                        UPDATE metadata_table SET meta_contents = %s WHERE meta_id = %s
                        """, (psycopg2.extras.Json(metadata_contents), meta_node_id))

                    self.db_connection.commit()

    def update_node_descendants_paths(self, new_parent_node_id):
        """
        This method updates node descendants paths
        :param new_parent_node_id: node_id of new parent node in tree_table which descendants we update
        """

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

        """.format(self.tree_table_name, new_parent_node_id)

        self.cursor.execute(recursive_update_command)
        self.db_connection.commit()

    def get_projections_on_path(self, path):
        """
        Method returns list of nodes on path
        :param path: path to node string
        :return: list of nodes paths strings
        """
        path = self.__split_path(path)

        assert isinstance(path, list), 'Path is not a list!'
        paths_request_command = """
        WITH node_on_path AS (
            SELECT node_id AS node_on_path_id FROM {0} WHERE path=%s::varchar[] AND projection_name=%s
        )
        SELECT {0}.name FROM {0}, node_on_path WHERE {0}.parent_id=node_on_path.node_on_path_id
        """.format(self.tree_table_name)

        self.cursor.execute(paths_request_command, (path, self.projection_id))
        return [row[0] for row in self.cursor]

    def move_projection(self, node_to_move_path, new_node_root_path):
        """
        This methods moves node in a tree to a new parent node, and updates node paths accordingly
        :param node_to_move_path: path to node which will be moved as list of strings
        :param new_node_root_path: path to new parent node list of strings
        """
        node_to_move_path = self.__split_path(node_to_move_path)
        new_node_root_path = self.__split_path(new_node_root_path)

        assert isinstance(node_to_move_path, list), 'Input old path is not a list!'
        assert isinstance(new_node_root_path, list), 'Input new path is not a list!'
        if node_to_move_path == new_node_root_path:
            raise Exception('Cannot move node onto itself!')

        node_ids_command = """
        WITH
        old_root_node AS (
            SELECT node_id AS old_root_node_id FROM {0} WHERE path=%s::varchar[]
        ),
        new_root_node AS (
            SELECT node_id AS new_root_node_id FROM {0} WHERE path=%s::varchar[]
        )
        SELECT * FROM old_root_node, new_root_node

        """.format(self.tree_table_name)

        self.cursor.execute(node_ids_command, (node_to_move_path, new_node_root_path))

        fetch_result = self.cursor.fetchone()

        if fetch_result is not None:
            node_id, new_parent_id = fetch_result
        else:
            raise OSError('Attempting to move non existant node on path: {}'.format(node_to_move_path))

        move_command = """
        UPDATE {0} SET parent_id={1} WHERE node_id = {2}
        """.format(self.tree_table_name, new_parent_id, node_id)
        self.cursor.execute(move_command)
        self.db_connection.commit()

        self.update_node_descendants_paths(new_parent_id)

    def remove_projection(self, node_path):
        """
        This method removes node from tree_table and all of it`s descendants
        :param node_path: path to node list of strings
        """
        node_path = self.__split_path(node_path)

        assert isinstance(node_path, list), 'Node path is not a list!'

        fetch_node_to_remove_id_command = """
        SELECT node_id FROM {0} WHERE path=%s::varchar[]
        """.format(self.tree_table_name)
        self.cursor.execute(fetch_node_to_remove_id_command, (node_path,))
        self.db_connection.commit()

        fetch_result = self.cursor.fetchone()

        if fetch_result is not None:
            node_to_remove_id = fetch_result[0]
        else:
            raise OSError('Attempting to delete non existant node on path {}'.format(node_path))

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
        """.format(self.tree_table_name, node_to_remove_id)

        self.cursor.execute(node_removal_command)
        self.db_connection.commit()

    def is_managing_path(self, path):
        """
        Check if projector is managing path
        :param path: projection path string
        :return: bool
        """
        # This command checks existance of projection row in tree_table by path
        # Run command and return check result as bool
        self.cursor.execute("""
        SELECT path FROM tree_table WHERE concat( '/', array_to_string(path[2:array_upper(path, 1)], '/'))=%s;
        """, (path,))

        return self.cursor.fetchone() is not None

    def get_attributes(self, path):
        """
        Get attributes of projection on given path
        :param path: path to projection string
        :return: attributes dict
        """
        attributes_order = ['st_atime', 'st_mtime', 'st_ctime', 'st_size', 'st_mode', 'st_nlink', 'st_ino']

        self.cursor.execute("""
        SELECT st_atime, st_mtime, st_ctime, st_size, st_mode, st_nlink, st_ino FROM get_projection_node_attributes(%s)
        """, (path,))
        # Fetching results of query
        attributes = self.cursor.fetchone()

        # Setting projection attributes dictionary using query results
        attributes = {el[0]: el[1] for el in zip(attributes_order, attributes)}

        # Setting appropriate types access modes for projection types for FUSE
        access_modes = {'file': (stat.S_IFREG | 0o0777),
                        'metadata': (stat.S_IFREG | 0o0777),
                        'directory': (stat.S_IFDIR | 0o0777)}

        attributes['st_mode'] = access_modes[attributes['st_mode']]
        return attributes

    def update_projection_size_attribute(self, path, size):
        path = self.__split_path(path)

        size_update_command = """
        WITH projection_on_path AS (
            SELECT node_id FROM {0} WHERE path = %s::varchar[]
        )
        UPDATE {1} SET st_size={2} FROM projection_on_path WHERE {1}.node_id=projection_on_path.node_id
        """.format(self.tree_table_name, self.projections_attributes_table_name, size)

        self.cursor.execute(size_update_command, (path,))
        self.db_connection.commit()

    def open_resource(self, path):
        """
        Opens resource on path and returns it`s header and contents stream
        :param path path string
        :return file_header
        :return resource_io
        """
        path = self.__split_path(path)

        projection_id_and_uri_query = """
        SELECT node_id, uri, type FROM {0} WHERE path = %s::varchar[]
        """.format(self.tree_table_name)

        self.cursor.execute(projection_id_and_uri_query, (path,))

        node_id, uri, projecton_type = self.cursor.fetchone()

        content = self.projection_driver.get_uri_contents_as_bytes(uri)
        logger.info('Got path content: %s\n', path)

        projection_size = len(content)

        update_projection_attributes_command = """
        UPDATE {0} SET st_size=%s WHERE node_id=%s
        """.format(self.projections_attributes_table_name)
        self.cursor.execute(update_projection_attributes_command, (projection_size, node_id,))

        file_header = 3
        resource_io = io.BytesIO(content)

        self.db_connection.commit()

        return file_header, resource_io