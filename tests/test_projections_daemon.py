#!/usr/bin/env python3

# Copyright 2016  Anton Bragin, Victor Svekolkin
#
# This file is part of Projections.
#
# Projections is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Projections is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Projections.  If not, see <http://www.gnu.org/licenses/>.

import logging
import logging.config
import os
import subprocess
import time
from unittest import TestCase

import psycopg2
import yaml

from projections_daemon import ProjectionsDaemon


class TestProjectionsDaemon(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_projections_logger')

        # Initializing database connection which will be used during tests
        with open('database_connection_config.yaml') as y_f:
            database_connection_parameters = yaml.safe_load(y_f)

        database_host = database_connection_parameters['database_host']
        database_port = database_connection_parameters['database_port']
        user_name = database_connection_parameters['user_name']
        user_password = database_connection_parameters['user_password']

        # Opening connection with database
        cls.db_connection = psycopg2.connect(database="projections_database",
                                             user=user_name,
                                             password=user_password,
                                             host=database_host,
                                             port=database_port)
        # Creating cursor, which will be used to interact with database
        cls.cursor = cls.db_connection.cursor()

        cls.cursor.execute(" DELETE FROM projections.projections WHERE projection_name~'test_projection'; ")
        cls.db_connection.commit()

    @classmethod
    def tearDownClass(cls):
        cls.cursor.execute(" DELETE FROM projections.projections WHERE projection_name~'test_projection'; ")
        cls.db_connection.commit()

        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        self.daemon = ProjectionsDaemon()

    def tearDown(self):
        # Stopping daemon and it`s projectors
        self.daemon.stop_daemon()

        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM projections.projections WHERE projection_name~'test_projection'; ")
        self.db_connection.commit()

    def test_get_projections(self):
        """
        Tests daemon get_projections returned list of projections correctness
        """
        for i in range(1, 6):
            self.cursor.execute("""
            INSERT INTO projections.projections (projection_name, mount_point, driver, prototype)
            VALUES ('test_projection_{0}', '/{0}', 'test_driver_{0}', 'test')
            RETURNING projection_id
            """.format(i))

            projection_id = self.cursor.fetchone()[0]

            self.db_connection.commit()

            reference_result = 'Projection id: {1} Projection name: test_projection_{0} ' \
                               'Mount point: /{0} Driver: test_driver_{0} Projector PID: None'.format(i,
                                                                                                      projection_id)
            self.assertIn(reference_result, self.daemon.get_projections(),
                          msg='Checking if daemon lists projections properly.')
            improper_result = 'Projection id: {1} Projection name: test_projection_{0} ' \
                              'Mount point: /{0} Driver: test_driver_{0} Projector PID: None'.format(i + 1,
                                                                                                     projection_id + 1)
            self.assertNotIn(improper_result, self.daemon.get_projections(),
                             msg='Checking if daemon not list misplaced projections.')

    def test_project(self):
        """
        Tests if daemon correctly initializes creation of projection
        """
        # Creating projection subprocess
        self.daemon.project('test_projection', 'tests/mnt', 'tests/projections_configs/test_metadata_operations.yaml',
                            'fs_driver')

        self.cursor.execute("""
        SELECT 1 FROM projections.projections WHERE projection_name = 'test_projection'
        """)

        self.assertIsNotNone(self.cursor.fetchone(), msg='Checking if projection where created by daemon.')

        self.cursor.execute("""
        SELECT projection_name, mount_point, driver FROM projections.projections WHERE projection_name = 'test_projection'
        """)

        projection_db_entry = self.cursor.fetchone()

        reference_entry = ['test_projection', 'tests/mnt', 'fs_driver']
        self.assertListEqual(reference_entry, list(projection_db_entry), msg='Checking projection entry correctness.')

    def test_stop(self):
        """
        Tests if daemon correctly stops underlying Projector subprocess
        """
        self.daemon.project('test_projection', 'tests/mnt', 'tests/projections_configs/test_metadata_operations.yaml',
                            'fs_driver')

        projectors = self.daemon.projections

        self.cursor.callproc('projections.daemon_get_projections')

        test_projection_id = None
        for row in self.cursor:
            # Dealing with psycopg2 inconsistency here, sometimes cursor returns ("",)
            if len(row) == 1:
                continue

            projection_id, projection_name, mount_point, driver, projector_pid, prototype = row
            if projection_name == 'test_projection':
                test_projection_id = projection_id

        self.assertTrue(any(projectors), msg='Checking if projectors dict is not empty after projection start')

        self.cursor.execute("""
        SELECT projector_pid FROM projections.projections WHERE projection_name = 'test_projection'
        """)
        projection_pid_before_stop = self.cursor.fetchone()[0]
        self.assertIsInstance(projection_pid_before_stop, int, msg='Checking if projection pid where setted.')

        # Stopping projector
        self.daemon.stop('test_projection')
        projectors = self.daemon.projections
        self.assertIsNone(projectors[test_projection_id]['projector_subprocess'],
                          msg='Checking if projection projector was set to None after stop.')

        self.cursor.execute("""
        SELECT projector_pid FROM projections.projections WHERE projection_name = 'test_projection'
        """)

        projection_pid_after_stop = self.cursor.fetchone()[0]
        self.assertIsNone(projection_pid_after_stop, msg='Checking if projector pid where set correctly.')

    def test_start(self):
        """
        Tests if daemon correctly starts underlying Projector subprocess
        """
        if os.path.ismount('tests/mnt'):
            raise OSError('Mount point tests/mnt is busy, please unmount it to contionue tests!')

        self.daemon.project('test_projection', 'tests/mnt', 'tests/projections_configs/test_metadata_operations.yaml',
                            'fs_driver')

        # Stopping projector
        self.daemon.stop('test_projection')

        projectors = self.daemon.projections

        self.cursor.callproc('projections.daemon_get_projections')

        self.logger.debug("Row")
        test_projection_id = None
        for row in self.cursor:
            # Dealing with psycopg2 inconsistency here, sometimes cursor returns ("",)
            if len(row) == 1:
                continue

            projection_id, projection_name, mount_point, driver, projector_pid, prototype = row
            self.logger.debug(projection_id)
            if projection_name == 'test_projection':
                test_projection_id = projection_id

        projector = projectors[test_projection_id]['projector_subprocess']

        self.assertIsNone(projector,
                          msg='Checking if projection projector was set to None after stop.')

        self.daemon.start('test_projection')

        projectors = self.daemon.projections

        projector = projectors[test_projection_id]['projector_subprocess']

        self.assertIsNotNone(projector,
                             msg='Checking if projection`s projector was set after start.')

        self.assertIsInstance(projector, subprocess.Popen, msg='Checking if projector is subprocess instance.')

    def test_remove_projection(self):
        """
        Tests if daemon correctly removes projection
        :return:
        """
        self.daemon.project('test_projection', 'tests/mnt', 'tests/projections_configs/test_metadata_operations.yaml',
                            'fs_driver')

        self.daemon.remove_projection('test_projection')

        # Checking if projection is in projections table
        self.cursor.execute("""
        SELECT 1 FROM projections.projections WHERE projection_name~'test_projection';
        """)
        is_projection_exists = self.cursor.fetchone() is None

        self.assertTrue(is_projection_exists,
                        msg='Checking if projection were removed properly.')

    def test_search(self):
        """
        Tests projection search by daemon
        """
        self.daemon.project('test_projection', 'tests/mnt', 'tests/projections_configs/test_metadata_operations.yaml',
                            'fs_driver')
        # Projector needs some time to initialize
        time.sleep(1.5)

        reference_result = {'/test_dir/fasta_file_1.fasta', '/test_dir/fasta_file_2.fasta',
                            '/test_dir/fasta_file_3.fasta', '/test_dir/fasta_file_4.fasta'}
        # Check search based on node properties
        search_result = self.daemon.search('test_projection', '/',
                                           query="""
                                            nodes.node_name ~ 'fasta$'
                                           """)
        self.logger.debug('Search results: %s', search_result)
        self.assertEqual(reference_result, set(search_result), msg='Checking if search returned proper results.')

        reference_result = {'/test_dir/bam_file_1.bam', '/test_dir/bam_file_2.bam', '/test_dir/bam_file_3.bam',
                            '/test_dir/bam_file_4.bam', '/test_dir/bam_file_5.bam'}
        # Check search based on metadata contents
        search_result = self.daemon.search('test_projection', '/',
                                           query="""
                                            nodes.node_name ~ 'bam$'
                                            AND nodes.node_id
                                                IN (SELECT head_node_id
                                                    FROM links
                                                    WHERE metadata_content->>'name'~'json$')
                                           """)
        self.logger.info('Search results: %s', search_result)
        self.assertEqual(reference_result, set(search_result), msg='Checking if search returned proper results.')

        for i in range(1, 6):
            # Check search on different path
            reference_result = {'/test_dir/sample_{}/vcf_file.vcf'.format(i),
                                '/test_dir/sample_{}/rerun/vcf_file.vcf'.format(i)}
            # Check search based on metadata contents
            search_result = self.daemon.search('test_projection', '/test_dir/sample_{}'.format(i),
                                               query="""
                                                nodes.node_name ~ 'vcf$'
                                               """)
            self.logger.debug('Search results: %s', search_result)
            self.assertEqual(reference_result, set(search_result), msg='Checking if search returned proper results.')
