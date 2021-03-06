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

import argparse
import copy
import logging
import logging.config
import os
import stat
import types

import objectpath
import psycopg2
import psycopg2.extras
import yaml

from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import PrototypeDeserializer

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('db_projector')


class DBProjector:
    """
    This class contains test interface to PostgresSQL database, which will be used to store and work with projections.
    Current implementation is subject of future rewrites and optimisations,
    """

    def __init__(self, projection_id, projection_driver, prototypes_tree, root_uri, restarting_projection=False):
        """
        This method initializes database which will hold projection tree and associated metadata
        :param projection_id: projection uri id integer
        :param projection_driver: projection driver instance
        :param prototypes_tree: prototypes tree instance
        :param root_uri: root uri of a projection
        """
        self.projection_id = projection_id
        self.prototypes_tree = prototypes_tree
        self.root_uri = root_uri

        # Initializing projection driver
        self.projection_driver = projection_driver

        with open('database_connection_config.yaml') as y_f:
            database_connection_parameters = yaml.safe_load(y_f)

        database_host = database_connection_parameters['database_host']
        database_port = database_connection_parameters['database_port']
        user_name = database_connection_parameters['user_name']
        user_password = database_connection_parameters['user_password']

        # Opening connection with database
        self.db_connection = psycopg2.connect(database="projections_database",
                                              user=user_name,
                                              password=user_password,
                                              host=database_host,
                                              port=database_port)

        # Setting connection mode of connection
        self.db_connection.autocommit = False

        # Creating cursor, which will be used to interact with database
        self.cursor = self.db_connection.cursor()

        if restarting_projection:
            logger.info('Restarting projection!')
        else:
            self.build_projection()
        logger.debug('Initialized projection: {}'.format(self.projection_id))

    def __del__(self):
        """
        This method closes cursor and database connection at call
        """
        self.cursor.close()
        self.db_connection.close()

    def build_projection(self):
        """
        This method builds projection
        """

        root_node_metadata_content = self.projection_driver.get_uri_contents_as_dict(self.root_uri)

        # Add root node to tables
        self.cursor.callproc("projections.add_projection_node", [self.projection_id,
                                                                 '',
                                                                 None,
                                                                 self.root_uri,
                                                                 'DIR',
                                                                 psycopg2.extras.Json(root_node_metadata_content),
                                                                 []])

        # After insertion cursor returns id of root node, which will be used in tree building
        new_parent_id = self.cursor.fetchone()

        # Setting root projection attributes
        self.cursor.callproc("projections.set_projection_node_attributes", [self.projection_id,
                                                                            new_parent_id,
                                                                            1,
                                                                            'DIR'])

        self.db_connection.commit()

        # Starting projections tree construction
        self.db_build_tree({'/': self.prototypes_tree}, root_projection_uri=self.root_uri, parent_id=new_parent_id)
        # Bind metadata to data
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
                content['content_uri'] = uri

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
                self.cursor.callproc("projections.add_projection_node", [self.projection_id,
                                                                         name,
                                                                         current_parent_id,
                                                                         uri,
                                                                         prototype_types_binding[prototype.type],
                                                                         psycopg2.extras.Json(node_metadata_content),
                                                                         prototype.meta_link])

                # Fetching inserted projection id which will be parent id for lower level projections
                new_parent_id = self.cursor.fetchone()

                # Setting projection attributes
                self.cursor.callproc("projections.set_projection_node_attributes", [self.projection_id, new_parent_id,
                                                                                    1,
                                                                                    prototype_types_binding[
                                                                                        prototype.type]])

                # If new_parent_id is None, when projection already exists in table
                if new_parent_id is not None:
                    self.db_build_tree(prototype.children, new_parent_id, root_projection_uri=uri)
                    # Commit all changes
                    self.db_connection.commit()

    def bind_metadata_to_data(self):
        """
        Performs data-metadata binding according to meta_links of nodes
        """
        try:
            self.cursor.callproc("projections.projector_bind_metadata_to_data", [self.projection_id])
            self.db_connection.commit()
        except psycopg2.ProgrammingError as e:
            logger.fatal('Error encountered during metadata-data binding: %s', e)

    def get_projections_on_path(self, path):
        """
        This method returns list of nodes on path
        :param path: path to node string
        :return: list of nodes paths strings
        """
        self.cursor.callproc("projections.get_projections_on_path", [self.projection_id, path])
        return self.cursor.fetchone()[0]

    def move_projection(self, node_to_move_path, new_node_root_path):
        """
        This methods moves node in a tree to a new parent node, and updates node paths accordingly
        :param node_to_move_path: path to node which will be moved as list of strings
        :param new_node_root_path: path to new parent node list of strings
        """

        self.cursor.callproc("projections.move_node_on_path", [self.projection_id,
                                                               node_to_move_path,
                                                               new_node_root_path])
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
        self.cursor.callproc("projections.remove_node_on_path", [self.projection_id, node_path])
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
        # This command checks existence of projection row in tree_table by path
        # Run command and return check result as bool
        self.cursor.callproc("projections.projector_is_managing_path", [self.projection_id, path])

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

    def open_resource(self, path):
        """
        Opens resource on path and returns it`s header and contents stream
        :param path path string
        :return file_header
        :return resource_io
        """
        self.cursor.callproc("projections.get_uri_of_projection_on_path", [self.projection_id, path])

        node_id, uri, node_on_path_st_mode = self.cursor.fetchone()

        content_context_manager = self.projection_driver.get_uri_contents_as_bytes(uri)
        logger.info('Got path content: %s\n', path)

        file_header = 3

        return file_header, content_context_manager

    def update_projection_size_attribute(self, path, projection_size):
        """

        :param path:
        :param projection_size:
        :return:
        """
        self.cursor.callproc("projections.get_uri_of_projection_on_path", [self.projection_id, path])
        node_id, _, node_on_path_st_mode = self.cursor.fetchone()

        # Updating projection attributes
        self.cursor.callproc("projections.set_projection_node_attributes", [self.projection_id,
                                                                            node_id,
                                                                            projection_size,
                                                                            node_on_path_st_mode])

        self.db_connection.commit()


def main(projection_name, projection_id, mount_point, prototype, driver, restart):
    # TODO consider driver addition from config file
    from drivers.aws_s3_driver import S3Driver
    from drivers.fs_driver import FSDriver
    from drivers.genbank_driver import GenbankDriver
    from drivers.iontorrent_driver import TorrentSuiteDriver
    from drivers.sra_driver import SRADriver

    drivers = {
        'iontorrent': TorrentSuiteDriver,
        'fs_driver': FSDriver,
        'genbank_driver': GenbankDriver,
        'sra_driver': SRADriver,
        'aws_s3_driver': S3Driver
    }

    mount_root, mount_point_name = os.path.split(mount_point)
    data_directory = os.path.join(mount_root, '.projections_contents', projection_name, mount_point_name)

    # Saving data to hidden directory near mount point
    if not os.path.exists(data_directory):
        os.makedirs(data_directory)

    projection_filesystem = ProjectionFilesystem(mount_point, data_directory)

    # Loading projection prototype and driver config
    projection_configuration = PrototypeDeserializer(prototype)

    projection_driver = drivers[driver](projection_configuration.resource_uri,
                                        projection_configuration.driver_config_path)
    # Initializing db projector
    projector = DBProjector(projection_id, projection_driver,
                            projection_configuration.prototype_tree,
                            projection_configuration.root_projection_uri,
                            restart)
    projection_filesystem.projection_manager = projector

    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    fuse = FUSE(projection_filesystem, mount_point, foreground=True, nonempty=True, nothreads=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DBProjector')

    parser.add_argument('-p_n', '--projection_name', required=True,
                        help='Id of a projection to be created (required).')

    parser.add_argument('-p_id', '--projection_id', required=True,
                        help='Id of a projection to be created (required).')
    parser.add_argument('-m', '--mount_point', required=True,
                        help='Folder in a system to mount projection to.')
    parser.add_argument('-p', '--prototype', required=True,
                        help='Path to prototype file to create projection for.')
    parser.add_argument('-d', '--driver', required=True,
                        help='Name of the driver to use with projection')
    parser.add_argument('-r', '--restart', action='store_true',
                        help='Restart projection')

    args = parser.parse_args()

    main(args.projection_name, args.projection_id, args.mount_point, args.prototype, args.driver, args.restart)
