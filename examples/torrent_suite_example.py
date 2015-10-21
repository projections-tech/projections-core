#!/usr/bin/env python3

import logging
import logging.config
import argparse
from fuse import FUSE
from projections import Projection, PrototypeDeserializer
from iontorrent import TorrentSuiteDriver, TorrentSuiteProjector
from filesystem import ProjectionFilesystem

def main(cfg_path, mountpoint, data_folder, foreground=True):
    LOGIN = 'ionadmin'
    PASSWORD = '0ECu1lW'

    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)

    # PrototypeDeserializer class builds prototype tree using which projections will be created,
    # specifies root projection uri and resource uri which will be projected
    projection_configuration = PrototypeDeserializer(cfg_path)
    # Root projection from which projection tree is build
    root_projection = Projection('/', projection_configuration.root_projection_uri)
    # This driver is used to connect with projection resource
    projection_driver = TorrentSuiteDriver(projection_configuration.resource_uri, LOGIN, PASSWORD)

    projection_filesystem.projection_manager = TorrentSuiteProjector(projection_driver,
                                                                     root_projection,
                                                                     projection_configuration.prototype_tree)

    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Torrent Suite projection example.')
    parser.add_argument('-m_p', '--mount_point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d_d', '--data_directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-cfg', '--config_path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory)
