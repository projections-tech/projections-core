#!/usr/bin/env python3

import argparse
import json
import io
import logging
import logging.config
import os
import time

from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import Projection, ProjectionDriver, ProjectionTree, Projector, PrototypeDeserializer

logger = logging.getLogger('filesystem_projection')


class FSDriver(ProjectionDriver):

    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """
        # Directory projection returns list of it`s children as metadata
        logger.debug(uri)
        if os.path.isdir(uri):
            return {'name': os.path.basename(uri), 'size': os.path.getsize(uri), 'type': 'dir',
                    'children': [os.path.join(os.path.abspath(uri), p) for p in os.listdir(uri)], 'extension': None}
        else:
            return {'name': os.path.basename(uri), 'size': os.path.getsize(uri), 'type': 'file',
                    'extension': os.path.splitext(uri)[1]}

    def get_uri_contents_as_stream(self, uri):
        """
        Load uri contents as stream
        :param uri: URI string
        :return: content bytes
        """
        if os.path.isfile(uri):
            with open(uri, 'rb') as f:
                return f.read()
        elif os.path.isdir(uri):
            return json.dumps(os.listdir(uri)).encode()


class FSProjector(Projector):
    def __init__(self, driver, root_projection, root_prototype):
        """
        Initializes filesystem Projector with driver, assigns root projection, builds prototype and projection tree.
        :param driver: instance of FSDriver
        :param prototype_tree: tree of ProjectionPrototype objects to build projection upon
        """
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree()
        self.root_projection = root_projection
        self.projection_tree.add_projection(self.root_projection, None)

        self.create_projection_tree({'/': root_prototype}, projection_tree=self.projection_tree, parent_projection=self.root_projection)
        self.projections = self.projection_tree.projections

    def is_managing_path(self, path):
        return path in self.projections

    def get_projections(self, path):
        logger.info('Requesting projections for path: %s', path)
        projections = []
        for p in self.projections:
            if p.startswith(path):
                # limit projections to one level only
                logging.debug('Analyzing projection candidate: %s', p)
                suffix = p[len(path):]
                if suffix and suffix[0] == '/':
                    suffix = suffix[1:]
                logger.debug('Path suffix: %s', suffix)
                if suffix and '/' not in suffix:
                    projections.append('/' + suffix)

        logger.debug('Returning projections: %s', projections)
        return projections

    def get_attributes(self, path):
        assert path in self.projections

        projection = self.projections[path]

        now = time.time()
        attributes = dict()

        # Set projection attributes

        # This is implementation specific and should be binded to projector data
        attributes['st_atime'] = now
        # This may be implemented as last projection cashing time is casing is enabled
        attributes['st_mtime'] = now
        # On Unix this is time for metedata modification we can use the same conception
        attributes['st_ctime'] = now
        # If this is projection the size is zero
        attributes['st_size'] = projection.size
        # Set type to link anf grant full access to everyone
        attributes['st_mode'] = (projection.type | 0o0777)
        # Set number of hard links to 0
        attributes['st_nlink'] = 0
        # Set id as inode number.
        attributes['st_ino'] = 1

        return attributes

    def open_resource(self, path):
        uri = self.projections[path].uri

        content = self.driver.get_uri_contents_as_stream(uri)
        logger.info('Got path content: %s\n', path)

        self.projections[path].size = len(content)

        file_header = 3
        resource_io = io.BytesIO(content)

        return file_header, resource_io


# For smoke testing
def main(cfg_path, mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)

    projection_configuration = PrototypeDeserializer(cfg_path)

    fs_driver = FSDriver()

    root_projection = Projection('/', projection_configuration.root_projection_uri)

    projection_filesystem.projection_manager = FSProjector(fs_driver, root_projection,
                                                           projection_configuration.prototype_tree)
    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    parser = argparse.ArgumentParser(description='Local filesystem projection.')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory)
