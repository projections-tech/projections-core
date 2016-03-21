import getpass
import logging
import logging.config
import os
import subprocess
import sys
import time
from unittest import TestCase

import psycopg2

from tests.mock import MockResource

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('genbank_test')

class TestGenbankProjectionContents(TestCase):
    @classmethod
    def setUpClass(cls):
        # Initializing database connection which will be used during tests
        cls.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))
        # Creating cursor, which will be used to interact with database
        cls.cursor = cls.db_connection.cursor()

    @classmethod
    def tearDownClass(cls):
        # Removing test projection entries from projections db
        cls.cursor.execute(" DELETE FROM projections_table WHERE projection_name='test_genbank_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        self.mock_resource = MockResource('tests/genbank_mock.json')

        subprocess.call([sys.executable,
                         'thin_daemon.py',
                         '--project',
                         '-p_t', 'genbank',
                         '-p_n', 'test_genbank_projection',
                         '-m', 'tests/mnt',
                         '-d', 'tests/data',
                         '-c', 'tests/test_genbank_config.yaml'])

        # Wait to initialize Projector
        time.sleep(0.7)

    def tearDown(self):
        subprocess.call([sys.executable,
                         'thin_daemon.py',
                         '-stop', 'test_genbank_projection'])

        # Unmounting projection dir
        subprocess.Popen(['fusermount', '-u', 'tests/mnt'])

        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM projections_table WHERE projection_name='test_genbank_projection' ")
        self.db_connection.commit()

        self.mock_resource.deactivate()

    def test_projections_contents(self):
        """
        Tests if file projections contents
        """

        logger.debug('Dir contents: %s', os.listdir('tests/mnt/'))

        # Loading test contents
        with open('tests/mock_resource/genbank_mock_data/mock_contents.txt') as f_f:
            test_contents = f_f.readlines()

        # Checking fasta file projection contents
        with open('tests/mnt/test_query_folder/sequence.fasta') as p_f:
            projected_fasta = []
            for line in p_f:
                logger.debug('Fasta line: %s', line)
                projected_fasta.append(line)

        logger.debug('Projected fasta length: %s', len(projected_fasta))

        self.assertEqual(projected_fasta, test_contents, msg='Check if fasta file contents loaded properly.')

        # Checking gb file projection contents
        with open('tests/mnt/test_query_folder/sequence.gb') as p_f:
            projected_gb = p_f.readlines()
        self.assertEqual(projected_gb, test_contents, msg='Check if gb file contents loaded properly.')
