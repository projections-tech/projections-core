#!/usr/bin/env python3

import argparse
import getpass
import logging
import logging.config
import os

import psycopg2

from db_projector import DBProjector
from drivers.iontorrent import TorrentSuiteDriver
from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import PrototypeDeserializer
from tests.mock import MockResource

logger = logging.getLogger('projection_daemon')

class ThinDaemon:
    def __init__(self, command_line_args):

        self.drivers = {
            'iontorrent': TorrentSuiteDriver
        }
        if command_line_args.project:
            logger.info('Daemon projecting!')

            self.projection_type = command_line_args.projection_type
            self.projection_config_path = command_line_args.config_path
            self.projection_mount_point = command_line_args.mount_point
            self.projection_data_directory = command_line_args.data_directory
            self.projection_name = command_line_args.projection_name

            self.run_projection()

        elif command_line_args.search:
            # Opening connection with database
            self.db_connection = psycopg2.connect(
                "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))

            self.cursor = self.db_connection.cursor()

            self.projection_name = command_line_args.projection_name

            self.perform_search(command_line_args.projection_name,
                                command_line_args.search_path,
                                command_line_args.search_query)

    def run_projection(self):
        # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
        # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/

        projection_filesystem = ProjectionFilesystem(self.projection_mount_point, self.projection_data_directory)

        projection_configuration = PrototypeDeserializer(self.projection_config_path)

        # For testing purposes use iontorrent config
        # TODO: implement driver configuration (PROJ-9)
        projection_driver = self.drivers[self.projection_type](projection_configuration.resource_uri,
                                                               'ionadmin', '0ECu1lW')

        projector = DBProjector(self.projection_name, projection_driver,
                                projection_configuration.prototype_tree,
                                projection_configuration.root_projection_uri)
        projection_filesystem.projection_manager = projector
        fuse = FUSE(projection_filesystem, self.projection_mount_point, foreground=True, nonempty=True)

    def stop_projection(self):
        raise NotImplementedError('Projection stopping is not implemented yet!')

    def delete_projection(self):
        raise NotImplementedError('Projection deletion is not implemented yet!')

    def list_projections(self):
        raise NotImplementedError('Projection listing is not implemented yet!')

    def perform_search(self, projection_name, path, query):
        """
        Perform search in projection using SQL as query language
        :param projection_name: name of projection on which to perform search
        :param path: path or level on which search is performed string
        :param query: SQL query which is used to filter projections
        :return: paths that adhere to search conditions into stdout
        """

        path = path.rstrip('/')

        if projection_name != 'global':
            self.cursor.execute("""
            SELECT node_id
            FROM tree_table
            WHERE (
                concat( '/', array_to_string(path[2:array_upper(path, 1)], '/')) = %s AND projection_name = %s
                )
            """, (path, projection_name))
        else:
            self.cursor.execute("""
            SELECT node_id
            FROM tree_table
            WHERE (
                concat( '/', array_to_string(path[2:array_upper(path, 1)], '/')) = %s
                )
            """, (path,))

        node_on_path_id = self.cursor.fetchone()[0]

        self.cursor.execute("""
        WITH descendants_table AS (
            WITH RECURSIVE

            descendants_ids AS (
                WITH RECURSIVE tree AS (
                    SELECT node_id, ARRAY[%(node_id)s]::integer[] AS ancestors
                    FROM tree_table WHERE parent_id = %(node_id)s

                    UNION ALL

                    SELECT tree_table.node_id, tree.ancestors || tree_table.parent_id
                    FROM tree_table, tree
                    WHERE tree_table.parent_id = tree.node_id
                )
                SELECT node_id FROM tree WHERE %(node_id)s = ANY(tree.ancestors)
            )

            SELECT tree_table.* FROM tree_table, descendants_ids WHERE tree_table.node_id = descendants_ids.node_id
        )
        SELECT concat( '/', array_to_string(descendants_table.path[2:array_upper(descendants_table.path, 1)], '/')) {query}
        """.format(query=query), {'node_id': node_on_path_id})
        # TODO add projections table to database and take mount point name from there
        for row in self.cursor:
            print('mount/' + row[0][1:])







def main():
    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    mock_torrent_suite = MockResource('tests/torrent_suite_mock.json')

    parser = argparse.ArgumentParser(description='Projection Daemon')

    parser.add_argument('-p', '--project', action='store_true', help='perform search')
    parser.add_argument('-s', '--search', action='store_true', help='perform projection')
    parser.add_argument('-p_t', '--projection-type', help='type of current projection')
    parser.add_argument('-p_n', '--projection-name', help='name of current projection')
    parser.add_argument('-m', '--mount-point', help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', help='specifies projection configuration YAML file path')
    parser.add_argument('-s_p', '--search-path', help='path on which to perform search')
    parser.add_argument('-s_q', '--search-query', help='search query')

    args = parser.parse_args()

    if args.project and args.projection_name is None and args.projection_type is None:
        parser.error("project action requires --projection-type and --projection-name.")
    elif args.search and args.projection_name is None:
        parser.error("search action requires --projection-name.")

    daemon = ThinDaemon(args)

if __name__ == '__main__':
    main()
