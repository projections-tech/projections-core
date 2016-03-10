#!/usr/bin/env python3

import argparse
import json
import logging
import logging.config
import os

from db_projector import DBProjector
from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import ProjectionDriver, PrototypeDeserializer

logger = logging.getLogger('filesystem_projection')


class FSDriver(ProjectionDriver):

    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """
        # Directory projection returns list of it`s children as metadata
        if os.path.isdir(uri):
            return {'name': os.path.basename(uri),
                    'resource_uri': os.path.abspath(uri),
                    'size': os.path.getsize(uri),
                    'type': 'dir',
                    'children': [
                        self.get_uri_contents_as_dict(os.path.join(os.path.abspath(uri), p)) for p in os.listdir(uri)
                        ],
                    'extension': None}
        else:
            return {'name': os.path.basename(uri),
                    'resource_uri': os.path.abspath(uri),
                    'size': os.path.getsize(uri),
                    'type': 'file',
                    'extension': os.path.splitext(uri)[1]}

    def get_uri_contents_as_bytes(self, uri):
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


# For smoke testing
def main(cfg_path, mountpoint, data_folder, projection_name, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)

    projection_configuration = PrototypeDeserializer(cfg_path)

    fs_driver = FSDriver()

    projection_filesystem.projection_manager = DBProjector(projection_name, fs_driver,
                                                           projection_configuration.prototype_tree,
                                                           projection_configuration.root_projection_uri)

    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    parser = argparse.ArgumentParser(description='Local filesystem projection.')

    parser.add_argument('-p', '--projection-name', required=True, help='name of current projection')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory, args.projection_name)
