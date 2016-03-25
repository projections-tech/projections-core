#!/usr/bin/env python3

import argparse
import logging
import logging.config
import os
import signal
import sys
import types

import objectpath
import psycopg2

from db_projector import DBProjector
# TODO consider driver addition from config file
from drivers.aws_s3_driver import S3Driver
from drivers.fs_driver import FSDriver
from drivers.genbank_driver import GenbankDriver
from drivers.iontorrent_driver import TorrentSuiteDriver
from drivers.sra_driver import SRADriver

from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import PrototypeDeserializer

from tests.mock import MockResource
import getpass


class ThinDaemon:
    def __init__(self, command_line_args, script_dir=None):

        self.logger = logging.getLogger('projection_daemon')

        self.daemon_pid = os.getpid()
        self.script_dir = script_dir

        self.drivers = {
            'iontorrent': TorrentSuiteDriver,
            'fs_projection': FSDriver,
            'genbank': GenbankDriver,
            'sra_projection': SRADriver,
            'aws_s3_driver': S3Driver
        }
        # Opening connection with database
        self.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))

        self.cursor = self.db_connection.cursor()

        self.cursor.execute("SELECT relname FROM pg_class WHERE relname = 'projections_table' ")
        projections_table_exists = self.cursor.fetchone()

        # If no projections table found, create required tables
        if not bool(projections_table_exists):
            self.cursor.execute("""
            CREATE TABLE projections_table (
                projection_name varchar PRIMARY KEY UNIQUE,
                mount_path varchar UNIQUE,
                projector_pid int
            );""")
            self.db_connection.commit()

        if command_line_args.project and script_dir is not None and os.path.dirname(os.path.realpath(__file__)) == '/':
            self.logger.info('Daemon projecting!')

            handlers = logging.root.handlers
            # Setting StreamHandler log level to critical in order to suppress console output
            handlers[0].setLevel(logging.CRITICAL)

            self.projection_type = command_line_args.projection_type
            self.projection_config_path = os.path.join(script_dir, command_line_args.config_path)
            self.projection_mount_point = os.path.join(script_dir, command_line_args.mount_point)
            self.projection_data_directory = os.path.join(script_dir, command_line_args.data_directory)
            self.projection_name = command_line_args.projection_name

            self.run_projection()

        elif command_line_args.search:
            self.projection_name = command_line_args.projection_name

            # self.perform_search(command_line_args.projection_name,
            #                     command_line_args.search_path,
            #                     command_line_args.search_query)

            self.perform_search_objectpath(command_line_args.projection_name,
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
        UPDATE projections_table SET projector_pid=Null, mount_path=Null WHERE projector_pid = %s
        """, (self.daemon_pid,))

        self.db_connection.commit()

        self.cursor.close()
        self.db_connection.close()

    def __split_path(self, path):
        """
        Split path string in parts preserving path root
        :param path: path string
        :return: list of path parts strings
        """
        temp = []
        while True:
            head, tail = os.path.split(path)
            if head == '/' or head == '':
                if tail == '':
                    return [head]
                else:
                    temp.append(tail)
                    if not head == '':
                        temp.append(head)
                    return temp[::-1]
            temp.append(tail)
            path = head

    def run_projection(self):
        """
        Runs projection
        """

        # Check if projection exists in projections_table
        self.cursor.execute("""
        SELECT EXISTS (SELECT 1 FROM projections_table WHERE projection_name=%s)
        """, (self.projection_name,))

        is_projection_exists = bool(self.cursor.fetchone()[0])

        # Check if projection is managed
        self.cursor.execute("""
        SELECT projector_pid FROM projections_table WHERE projection_name=%s
        """, (self.projection_name,))
        is_projection_managed = self.cursor.fetchone()

        if is_projection_managed is not None:
            if bool(is_projection_managed[0]):
                sys.exit('Error: projection "{}" is already running!'.format(self.projection_name))

        if not is_projection_exists:
            # Add projection if it doesnt exists
            self.cursor.execute("""
            INSERT INTO projections_table (projection_name, mount_path, projector_pid)
            VALUES (%(pr_name)s, %(mount_pt)s, %(pr_pid)s)
            """, {'pr_name': self.projection_name, 'mount_pt': self.projection_mount_point, 'pr_pid': self.daemon_pid})
        else:
            # If projection is in projections_database update it`s projector_pid
            self.cursor.execute("""
            UPDATE projections_table SET projector_pid=%s, mount_path=%s WHERE projection_name = %s
            """, (self.daemon_pid, self.projection_mount_point, self.projection_name))

        self.db_connection.commit()

        # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
        # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/

        projection_filesystem = ProjectionFilesystem(self.projection_mount_point, self.projection_data_directory)

        # Loading projection prototype and driver config
        projection_configuration = PrototypeDeserializer(self.projection_config_path)

        projection_driver = self.drivers[self.projection_type](projection_configuration.resource_uri,
                                                               projection_configuration.driver_config_path,
                                                               self.script_dir)
        # Initializing db projector
        projector = DBProjector(self.projection_name, projection_driver,
                                projection_configuration.prototype_tree,
                                projection_configuration.root_projection_uri)
        projection_filesystem.projection_manager = projector

        fuse = FUSE(projection_filesystem, self.projection_mount_point, foreground=True, nonempty=True, nothreads=True)

    def stop_projection(self, projection_name):
        """
        Stops projection process
        :param projection_name: name of projection to stop string
        """
        # Fetch projector pid by projection name
        self.cursor.execute("""
        SELECT projector_pid FROM projections_table WHERE projection_name = %s
        """, (projection_name,))
        projector_pid = self.cursor.fetchone()
        if projector_pid is not None:
            if projector_pid[0] is None:
                self.logger.info('Projection "%s" is not running!', projection_name)
            else:
                self.logger.info('Stopping projection "%s"!', projection_name)
                self.logger.debug('Projector pid:%s', projector_pid)
                os.kill(projector_pid[0], signal.SIGTERM)
        else:
            self.logger.info('Projection "{}" does not exist!'.format(projection_name))

    def delete_projection(self, projection_name):
        """
        This method performs deletion of projection
        :param projection_name:
        """
        self.cursor.execute("""
        SELECT projector_pid FROM projections_table WHERE projection_name = %s
        """, (projection_name,))

        is_projection_running = bool(self.cursor.fetchone()[0])

        if is_projection_running:
            self.logger.info('Attempting do delete running projection!')
            status = input('Continue? y/n ')
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
        self.logger.info('Removing projection "%s"!', projection_name)

    def list_projections(self):
        """
        Lists projections in projections database
        """
        self.cursor.execute("""
        SELECT * FROM projections_table
        """)

        for row in self.cursor:
            self.logger.info('Projection name: {0}\tMount point: {1}\tProjector pid: {2}'.format(*row))

    def perform_search(self, projection_name, path, query):
        """
        Perform search in projection using SQL as query language
        :param projection_name: name of projection on which to perform search
        :param path: path or level on which search is performed string
        :param query: SQL query which is used to filter projections
        :return: paths that adhere to search conditions into stdout
        """
        if not path == '/':
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

    def perform_search_objectpath(self, projection_name, path, query):
        """
        Perform search in projection using ObjectPath as query language
        :param projection_name: name of projection on which to perform search
        :param path: path or level on which search is performed string
        :param query: ObjectPath query which is used to filter projections
        :return: paths that adhere to search conditions into stdout
        """

        if not path == '/':
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
        # Setting root node on which to perform search
        node_on_path_id = self.cursor.fetchone()[0]

        node_children_metadata = {'projections': []}

        def recursive_descendants_query(node_id, json_dict):
            self.cursor.execute(" SELECT ARRAY(SELECT node_id FROM tree_table WHERE parent_id=%s) ", (node_id,))
            children = self.cursor.fetchone()[0]

            self.cursor.execute(" SELECT name FROM tree_table WHERE node_id=%s", (node_id,))

            node_name = self.cursor.fetchone()[0]
            if node_name == '/':
                node_name = 'root'

            self.cursor.execute("""
            WITH node_metadata AS (
                SELECT node_id, meta_contents FROM metadata_table WHERE parent_node_id=%s
                )
            SELECT tree_table.name, node_metadata.meta_contents
            FROM tree_table, node_metadata
            WHERE node_metadata.node_id = tree_table.node_id
            """, (node_id,))

            node_metadata = dict()
            for row in self.cursor:
                node_metadata = row[1]

            json_dict['projections'].append({'name': node_name, 'node_id': node_id, 'metadata': node_metadata})
            for child in children:
                recursive_descendants_query(child, json_dict)

        recursive_descendants_query(node_on_path_id, node_children_metadata)

        tree = objectpath.Tree(node_children_metadata)
        # Query example: $.projections[@.metadata.status is "complete"]
        query = tree.execute(query)

        # Object path sometimes returns generator if user uses selectors, for consistency expand it using
        # list comprehension
        if isinstance(query, types.GeneratorType):
            query = [el for el in query]

        res_nodes_ids = [el['node_id'] for el in query]

        self.cursor.execute("""
        SELECT mount_path FROM projections_table WHERE projection_name=%s
        """, (self.projection_name,))

        mount_path = self.cursor.fetchone()[0]

        self.cursor.execute("""
        SELECT concat( '/', array_to_string(path[2:array_upper(path, 1)], '/'))
        FROM tree_table
        WHERE node_id = ANY(%s)
        """, (res_nodes_ids,))

        for row in self.cursor:
            print(os.path.join(mount_path, row[0].lstrip('/')))

    def bind_metadata_to_path(self, target_path, metadata_path):
        """
        This method binds metadata object on path to target path
        :param target_path: path to target with which metadata is binded, list of stings
        :param metadata_path: path to metadata object which will be binded, list of strings
        """

        target_path = self.__split_path(target_path)
        metadata_path = self.__split_path(metadata_path)

        binding_command = """
        WITH
        target AS (
            SELECT node_id AS target_id FROM tree_table WHERE path=%s::varchar[]
        ),
        metadata AS (
            SELECT node_id AS metadata_id FROM tree_table WHERE path=%s::varchar[]
        )
        INSERT INTO metadata_table (parent_node_id, node_id) SELECT target.target_id, metadata.metadata_id FROM target, metadata
        WHERE NOT EXISTS (SELECT * FROM target, metadata, metadata_table
        WHERE metadata_table.node_id = metadata.metadata_id AND metadata_table.parent_node_id = target.target_id)
        """
        self.cursor.execute(binding_command, (target_path, metadata_path))

        self.db_connection.commit()


def main():
    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

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

    mock_resource = MockResource('tests/torrent_suite_mock.json')

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
        os.chdir('/')
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
