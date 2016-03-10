#!/usr/bin/env python3

import argparse

from filesystem import ProjectionFilesystem
from fs_projection import FSDriver
from fuse import FUSE
from projections import Projector, PrototypeDeserializer


def main(cfg_path, mountpoint, data_folder, foreground=True):
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)

    # PrototypeDeserializer class builds prototype tree, using which projections will be created,
    # specifies root projection uri and resource uri which will be projected
    projection_configuration = PrototypeDeserializer(cfg_path)

    # This driver is used to connect with projection resource
    projection_driver = FSDriver()

    projection_filesystem.projection_manager = Projector(projection_driver,
                                                         projection_configuration.root_projection_uri,
                                                         projection_configuration.prototype_tree).projection_tree

    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filesystem projection example.')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory)

