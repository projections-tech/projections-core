import argparse
import logging
import logging.config
import os

from db_projector import DBProjector
from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import PrototypeDeserializer
from tests.mock import MockResource


class ThinDaemon:
    def __init__(self, action, **kwargs):
        pass

    def run_projection(self):
        # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
        # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
        projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)

        mock_torrent_suite = MockResource('tests/torrent_suite_mock.json')
        projection_configuration = PrototypeDeserializer(cfg_path)
        projection_driver = TorrentSuiteDriver(projection_configuration.resource_uri, 'ionadmin', '0ECu1lW')

        ion_torrent_projection_tree = DBProjector(projection_name, projection_driver,
                                                  projection_configuration.prototype_tree,
                                                  projection_configuration.root_projection_uri)
        projection_filesystem.projection_manager = ion_torrent_projection_tree
        fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
        return fuse


def main(cfg_path, mountpoint, data_folder, projection_name, foreground=True):
    d = ThinDaemon


if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    parser = argparse.ArgumentParser(description='Torrent Suite projection.')
    parser.add_argument('-p', '--projection-name', required=True, help='name of current projection')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory, args.projection_name)
