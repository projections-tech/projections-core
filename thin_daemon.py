#!/usr/bin/env python3

import argparse
import getpass
import logging
import logging.config
import os
import signal
import sys

import psycopg2

from db_projector import DBProjector
from drivers.iontorrent import TorrentSuiteDriver
from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import PrototypeDeserializer
from tests.mock import MockResource

logger = logging.getLogger('projection_daemon')

class ThinDaemon:
    def __init__(self, command_line_args, script_dir=None):

        self.daemon_pid = os.getpid()

        self.drivers = {
            'iontorrent': TorrentSuiteDriver
        }
        # Opening connection with database
        self.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))

        self.cursor = self.db_connection.cursor()

        self.cursor.execute("SELECT relname FROM pg_class WHERE relname = 'projections_table' ")
        is_projections_table_exist = self.cursor.fetchone()

        # If no projections table found, create required tables
        if not bool(is_projections_table_exist):
            self.cursor.execute("""
            CREATE TABLE projections_table (
                projection_name varchar PRIMARY KEY UNIQUE,
                mount_path varchar UNIQUE,
                projector_pid int
            );""")
            self.db_connection.commit()

        if command_line_args.project and script_dir is not None:
            logger.info('Daemon projecting!')
            logging.disable(logging.CRITICAL)
            self.projection_type = command_line_args.projection_type
            self.projection_config_path = os.path.join(script_dir, command_line_args.config_path)
            self.projection_mount_point = os.path.join(script_dir, command_line_args.mount_point)
            self.projection_data_directory = os.path.join(script_dir, command_line_args.data_directory)
            self.projection_name = command_line_args.projection_name

            # Opening connection with database
            self.db_connection = psycopg2.connect(
                "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))

            self.cursor = self.db_connection.cursor()

            self.run_projection()

        elif command_line_args.search:
            self.projection_name = command_line_args.projection_name

            self.perform_search(command_line_args.projection_name,
                                command_line_args.search_path,
                                command_line_args.search_query)

        elif command_line_args.stop_projection:
            self.stop_projection(command_line_args.stop_projection)

        elif command_line_args.delete_projection:
            self.delete_projection(command_line_args.delete_projection)

        elif command_line_args.list_projections:
            self.list_projections()

    def __del__(self):
        # Remove current pid from projections table on deconstruction
        self.cursor.execute("""
        UPDATE projections_table SET projector_pid=Null WHERE projector_pid = %s
        """, (self.daemon_pid,))

        self.db_connection.commit()

        self.cursor.close()
        self.db_connection.close()

    def run_projection(self):
        """
        Runs projection
        """

        # Check if projection exists in projections_table
        self.cursor.execute("""
        SELECT EXISTS (SELECT 1 FROM projections_table WHERE projection_name=%s)
        """, (self.projection_name,))

        is_projection_exists = bool(self.cursor.fetchone()[0])

        if not is_projection_exists:
            # Add projection if it doesnt exists
            self.cursor.execute("""
            INSERT INTO projections_table (projection_name, mount_path, projector_pid)
            VALUES (%(pr_name)s, %(mount_pt)s, %(pr_pid)s)
            """, {'pr_name': self.projection_name, 'mount_pt': self.projection_mount_point, 'pr_pid': self.daemon_pid})
        else:
            # If projection is in projections_database update it`s projector_pid
            self.cursor.execute("""
            UPDATE projections_table SET projector_pid=%s WHERE projection_name = %s
            """, (self.daemon_pid, self.projection_name))

        self.db_connection.commit()

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

    def stop_projection(self, projection_name):
        """
        Stops projection process
        :param projection_name: name of projection to stop string
        """
        # Fetch projector pid by projection name
        self.cursor.execute("""
        SELECT projector_pid FROM projections_table WHERE projection_name = %s
        """, (projection_name,))
        projector_pid = self.cursor.fetchone()[0]
        if projector_pid is None:
            logger.info('Projection "%s" is not running!', projection_name)
        else:
            logger.info('Stopping projection "%s"!', projection_name)
            os.kill(projector_pid, signal.SIGTERM)

    def delete_projection(self, projection_name):
        """
        This method performs deletion of projection
        :param projection_name:
        :return:
        """
        self.cursor.execute("""
        SELECT projector_pid FROM projections_table WHERE projection_name = %s
        """, (projection_name,))

        is_projection_running = bool(self.cursor.fetchone()[0])

        if is_projection_running:
            logger.info('Attempting do delete running projection!')
            status = input('Continue? y/n')
            if status == 'y':
                self.stop_projection(projection_name)
                self.cursor.execute("""
                DELETE FROM projections_table WHERE projection_name = %s
                """, (projection_name,))
                self.db_connection.commit()
            else:
                return None
        else:
            self.cursor.execute("""
            DELETE FROM projections_table WHERE projection_name = %s
            """, (projection_name,))
            self.db_connection.commit()
        logger.info('Removing projection "%s"!', projection_name)

    def list_projections(self):
        """
        Lists projections in projections database
        :return:
        """
        self.cursor.execute("""
        SELECT * FROM projections_table
        """)

        for row in self.cursor:
            logger.info('Projection name: {0}\tMount point: {1}\t'.format(*row[:2]))

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
        SELECT mount_path FROM projections_table WHERE projection_name=%s
        """, (self.projection_name,))

        mount_path = self.cursor.fetchone()[0]

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

        for row in self.cursor:
            print(os.path.join(mount_path, row[0].lstrip('/')))


def main():
    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    mock_torrent_suite = MockResource('tests/torrent_suite_mock.json')

    parser = argparse.ArgumentParser(description='Projection Daemon')

    parser.add_argument('-p', '--project', action='store_true', help='perform search')
    parser.add_argument('-s', '--search', action='store_true', help='perform projection')
    parser.add_argument('-stop', '--stop_projection', help='stop projection')
    parser.add_argument('-delete', '--delete_projection', help='delete projection')
    parser.add_argument('-list', '--list_projections', action='store_true', help='list projections')

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

    # Daemon creation code
    if args.project:
        try:
            pid = os.fork()
            if pid > 0:
                # exit first parent
                sys.exit(0)
        except OSError as e:
            sys.exit(1)

        # decouple from parent environment
        os.chdir("/")
        os.setsid()
        os.umask(0)

        # do second fork
        try:
            pid = os.fork()
            if pid > 0:
                # exit from second parent
                sys.exit(0)
        except OSError as e:
            sys.exit(1)

        daemon = ThinDaemon(args, script_dir)
    daemon = ThinDaemon(args)

if __name__ == '__main__':
    main()
