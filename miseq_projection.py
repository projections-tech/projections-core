#!/usr/bin/env python3

import argparse
import json
import logging
import logging.config
import os
import subprocess
import shutil
import time

from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import ProjectionDriver, Projector, PrototypeDeserializer

logger = logging.getLogger('miseq_projection')


class MiSeqDriver(ProjectionDriver):

    def __init__(self, username, shared_dir_path):
        self.shared_dir_path = shared_dir_path
        # Naming temp directory using seconds since epoch to get unique name
        self.temp_dir_path = ''.join(['.', 'temp_', str((time.time())).replace('.', '_')])
        # Creating temp directory into which shared directory will be mounted
        os.mkdir(self.temp_dir_path)
        # Mounting shared directory into temp dir
        subprocess.call(['sudo', 'mount', '-t',
                         'cifs', '-o', 'username={0}'.format(username),
                         shared_dir_path, self.temp_dir_path])

    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """
        # Return contents of temporary directory if uri is path to shared directory, else join path to temp dir and uri
        if uri == self.shared_dir_path:
            return {'name': os.path.basename(self.temp_dir_path),
                    'size': os.path.getsize(self.temp_dir_path), 'type': 'dir',
                    'children': [os.path.join(os.path.abspath(self.temp_dir_path), p) for p in os.listdir(self.temp_dir_path)],
                    'extension': None}
        else:
            uri = os.path.join(self.temp_dir_path, uri)

        # Directory projection returns list of it`s children as metadata
        if os.path.isdir(uri):
            return {'name': os.path.basename(uri), 'size': os.path.getsize(uri), 'type': 'dir',
                    'children': [os.path.join(os.path.abspath(uri), p) for p in os.listdir(uri)], 'extension': None}
        else:
            return {'name': os.path.basename(uri), 'size': os.path.getsize(uri), 'type': 'file',
                    'extension': os.path.splitext(uri)[1]}

    def get_uri_contents_as_bytes(self, uri):
        """
        Load uri contents as bytes massive
        :param uri: URI string
        :return: content bytes
        """
        if os.path.isfile(uri):
            with open(uri, 'rb') as f:
                return f.read()
        elif os.path.isdir(uri):
            return json.dumps(os.listdir(uri)).encode()

    def __del__(self):
        # Unmounting temp directory
        subprocess.call(['sudo', 'umount', '-l', self.temp_dir_path])
        # Removing temp directory after it is unmounted
        shutil.rmtree(self.temp_dir_path, ignore_errors=True)


# For smoke testing
def main(cfg_path, mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)

    projection_configuration = PrototypeDeserializer(cfg_path)

    miseq_driver = MiSeqDriver('vsvekolkin', projection_configuration.resource_uri)

    projection_filesystem.projection_manager = Projector(miseq_driver, projection_configuration.root_projection_uri,
                                                         projection_configuration.prototype_tree).projection_tree
    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    parser = argparse.ArgumentParser(description='MiSeq projection.')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory)

