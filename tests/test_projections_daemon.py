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

import getpass
import logging
import logging.config
import subprocess
import time
from unittest import TestCase

import psycopg2

from projections_daemon import ProjectionsDaemon


class TestProjectionsDaemon(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_projections_logger')

        # Initializing database connection which will be used during tests
        cls.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))
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

        projections_ids = []
        for i in range(1, 6):
            self.cursor.execute("""
            INSERT INTO projections.projections (projection_name, mount_point, driver)
            VALUES ('test_projection_{0}', '/{0}', 'test_driver_{0}')
            RETURNING projection_id
            """.format(i))
            projections_ids.append(self.cursor.fetchone()[0])

            self.db_connection.commit()

        reference_result = []
        for i, proj_id in enumerate(projections_ids):
            reference_result.append('Projection id: {1} Projection name: test_projection_{0} '
                                    'Mount point: /{0} Driver: test_driver_{0} Projector PID: None'.format(i + 1,
                                                                                                           proj_id))
        self.assertListEqual(reference_result, self.daemon.get_projections(),
                             msg='Checking id daemon lists projections properly.')

    def test_project(self):
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
        self.daemon.project('test_projection', 'tests/mnt', 'tests/projections_configs/test_metadata_operations.yaml',
                            'fs_driver')

        projectors = self.daemon.projectors
        projection_id = list(projectors.keys())[0]

        self.assertTrue(any(projectors), msg='Checking if projectors dict is not empty after projection start')

        self.cursor.execute("""
        SELECT projector_pid FROM projections.projections WHERE projection_name = 'test_projection'
        """)
        projection_pid_before_stop = self.cursor.fetchone()[0]
        self.assertIsInstance(projection_pid_before_stop, int, msg='Checking if projection pid where setted.')

        # Stopping projector
        self.daemon.stop('test_projection')
        projectors = self.daemon.projectors
        self.assertIsNone(projectors[projection_id]['projector_subprocess'],
                          msg='Checking if projection projector was set to None after stop.')

        self.cursor.execute("""
        SELECT projector_pid FROM projections.projections WHERE projection_name = 'test_projection'
        """)

        projection_pid_after_stop = self.cursor.fetchone()[0]
        self.assertIsNone(projection_pid_after_stop, msg='Checking if projector pid where set correctly.')

    def test_start(self):
        self.daemon.project('test_projection', 'tests/mnt', 'tests/projections_configs/test_metadata_operations.yaml',
                            'fs_driver')

        # Stopping projector
        self.daemon.stop('test_projection')

        projectors = self.daemon.projectors
        projection_id = list(projectors.keys())[0]

        projector = projectors[projection_id]['projector_subprocess']

        self.assertIsNone(projector,
                          msg='Checking if projection projector was set to None after stop.')

        self.daemon.start('test_projection')

        projectors = self.daemon.projectors
        projection_id = list(projectors.keys())[0]

        projector = projectors[projection_id]['projector_subprocess']

        self.assertIsNotNone(projector,
                             msg='Checking if projection`s projector was set after start.')

        self.assertIsInstance(projector, subprocess.Popen, msg='Checking if projector is subprocess instance.')

    def test_remove_projection(self):
        self.daemon.project('test_projection', 'tests/mnt', 'tests/projections_configs/test_metadata_operations.yaml',
                            'fs_driver')

        self.daemon.remove_projection('test_projection')

        # Checking if projection is in projections table
        self.cursor.execute("""
        SELECT 1 FROM projections.projections WHERE projection_name~'test_projection';
        """)
        is_projection_exists = self.cursor.fetchone()[0] is not None

        self.assertTrue(is_projection_exists,
                        msg='Checking if projection were removed properly.')

    def test_search(self):
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

        reference_result = {'/test_dir/bam_file_1.bam', '/test_dir/bam_file_2.bam', '/test_dir/bam_file_3.bam'}
        # Check search based on metadata contents
        search_result = self.daemon.search('test_projection', '/',
                                           query="""
                                            nodes.node_name ~ 'bam$'
                                            AND nodes.node_id
                                                IN (SELECT head_node_id
                                                    FROM links
                                                    WHERE metadata_content->'size'='23')
                                           """)
        self.logger.debug('Search results: %s', search_result)
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
