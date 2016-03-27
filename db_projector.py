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

    def build_projection(self):
        """
        This method builds projection
        """
        root_node_metadata_content = self.projection_driver.get_uri_contents_as_dict(self.root_uri)

        # Add root node to tables
        self.cursor.execute("""
        SELECT projections.add_projection_node(
            %(projection_id)s,
            %(node_name)s,
            %(parent_id)s,
            %(uri)s,
            %(type)s,
            %(metadata_content)s,
            %(meta_links)s
            )
        """, {'projection_id': self.projection_id,
              'parent_id': None,
              'node_name': '',
              'type': 'DIR',
              'uri': self.root_uri,
              'metadata_content': psycopg2.extras.Json(root_node_metadata_content),
              'meta_links': []})

        # After insertion cursor returns id of root node, which will be used in tree building
        new_parent_id = self.cursor.fetchone()

        # Setting root projection attributes
        self.cursor.execute("""
        SELECT projections.set_projection_node_attributes(%(tree_id)s, %(node_id)s, %(node_size)s, %(node_mode)s )
        """, {'tree_id': self.projection_id,
              'node_id': new_parent_id,
              'node_size': 1,
              'node_mode': 'DIR'})

        self.db_connection.commit()

        # Starting projections tree construction
        self.db_build_tree({'/': self.prototypes_tree}, root_projection_uri=self.root_uri, parent_id=new_parent_id)
        self.bind_metadata_to_data()

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

                node_metadata_content = self.projection_driver.get_uri_contents_as_dict(uri)

                # Inserting node into projections tree on completion this command returns inserted node id
                # which will be passed to lower level nodes
                self.cursor.execute("""
                SELECT projections.add_projection_node(
                        %(projection_id)s, %(node_name)s, %(parent_id)s,
                        %(uri)s, %(type)s, %(metadata_content)s, %(meta_links)s
                        )
                """, {'projection_id': self.projection_id,
                      'parent_id': current_parent_id,
                      'node_name': name,
                      'type': prototype_types_binding[prototype.type],
                      'uri': self.root_uri,
                      'metadata_content': psycopg2.extras.Json(node_metadata_content),
                      'meta_links': prototype.meta_link
                      })

                # Fetching inserted projection id which will be parent id for lower level projections
                new_parent_id = self.cursor.fetchone()

                # Setting projection attributes
                self.cursor.execute("""
                SELECT projections.set_projection_node_attributes(%(tree_id)s, %(node_id)s, %(node_size)s, %(node_mode)s )
                """, {'tree_id': self.projection_id,
                      'node_id': new_parent_id,
                      'node_size': 1,
                      'node_mode': prototype_types_binding[prototype.type]}
                                    )

                # If new_parent_id is None, when projection already exists in table
                if new_parent_id is not None:
                    self.db_build_tree(prototype.children, new_parent_id, root_projection_uri=uri)
                    # Commit all changes
                    self.db_connection.commit()

    def bind_metadata_to_data(self):
        """
        Performs data-metadata binding according to meta_link
        """
        self.cursor.execute("""
        SELECT projections.projector_bind_metadata_to_data(%s)
        """, (self.projection_id,))
        self.db_connection.commit()

    def get_projections_on_path(self, path):
        """
        Method returns list of nodes on path
        :param path: path to node string
        :return: list of nodes paths strings
        """
        self.cursor.execute("""
        SELECT projections.get_projections_on_path(%s, %s)
        """, (self.projection_id, path))
        return self.cursor.fetchone()[0]

    def move_projection(self, node_to_move_path, new_node_root_path):
        """
        This methods moves node in a tree to a new parent node, and updates node paths accordingly
        :param node_to_move_path: path to node which will be moved as list of strings
        :param new_node_root_path: path to new parent node list of strings
        """

        self.cursor.execute("""
        SELECT projections.move_node_on_path(%(tree_id)s, %(node_to_move_path)s, %(new_node_root_path)s )
        """, {'tree_id': self.projection_id,
              'node_to_move_path': node_to_move_path,
              'new_node_root_path': new_node_root_path}
                            )
        is_node_moved = self.cursor.fetchone()[0]
        if is_node_moved:
            self.db_connection.commit()
        else:
            raise OSError('Node on path:%s move to path:%s attempt failed!', node_to_move_path, new_node_root_path)

    def remove_projection(self, node_path):
        """
        This method removes node from tree_table and all of it`s descendants
        :param node_path: path to node list of strings
        """
        self.cursor.execute("""
        SELECT projections.remove_node_on_path(%s, %s)
        """, (self.projection_id, node_path))
        # Checking if node removal where successful.
        is_successful = self.cursor.fetchone()[0]
        if is_successful:
            logger.info('Succesfully removed projection on path: %s', node_path)
            self.db_connection.commit()
        else:
            raise OSError('Attempting to remove non existant projection on path: %s', node_path)

    def is_managing_path(self, path):
        """
        Check if projector is managing path
        :param path: projection path string
        :return: bool
        """
        # This command checks existance of projection row in tree_table by path
        # Run command and return check result as bool
        self.cursor.execute("""
        SELECT projections.projector_is_managing_path(%s, %s)
        """, (self.projection_id, path))

        return self.cursor.fetchone()[0]

    def get_attributes(self, path):
        """
        Get attributes of projection on given path
        :param path: path to projection string
        :return: attributes dict
        """
        attributes_order = ['st_atime', 'st_mtime', 'st_ctime', 'st_size', 'st_mode', 'st_nlink', 'st_ino']

        self.cursor.execute("""
        SELECT st_atime, st_mtime, st_ctime, st_size, st_mode, st_nlink, st_ino
        FROM projections.get_projection_node_attributes(%s, %s)
        """, (self.projection_id, path))
        # Fetching results of query
        attributes = self.cursor.fetchone()

        # Setting projection attributes dictionary using query results
        attributes = {el[0]: el[1] for el in zip(attributes_order, attributes)}

        # Setting appropriate types access modes for projection types for FUSE
        access_modes = {'REG': (stat.S_IFREG | 0o0777),
                        'DIR': (stat.S_IFDIR | 0o0777)}

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
        self.cursor.execute("""
        SELECT node_on_path_id, node_on_path_uri, node_on_path_st_mode
        FROM projections.get_uri_of_projection_on_path(%s, %s)
        """, (self.projection_id, path))

        node_id, uri, node_on_path_st_mode = self.cursor.fetchone()

        content = self.projection_driver.get_uri_contents_as_bytes(uri)
        logger.info('Got path content: %s\n', path)

        projection_size = len(content)

        # Updating projection attributes
        self.cursor.execute("""
        SELECT projections.set_projection_node_attributes(%(tree_id)s, %(node_id)s, %(node_size)s, %(node_mode)s )
        """, {'tree_id': self.projection_id,
              'node_id': node_id,
              'node_size': projection_size,
              'node_mode': node_on_path_st_mode})

        file_header = 3
        resource_io = io.BytesIO(content)

        self.db_connection.commit()

        return file_header, resource_io