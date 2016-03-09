import getpass
import logging
import logging.config
import os
import subprocess
import sys
import time
from unittest import TestCase

import psycopg2

import genbank
from db_projector import DBProjector
from projections import PrototypeDeserializer
from tests.mock import MockResource

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('genbank_test')


class TestGenbankProjector(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mock_resource = MockResource('tests/genbank_mock.json')

        # Initializing database connection which will be used during tests
        cls.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))
        # Creating cursor, which will be used to interact with database
        cls.cursor = cls.db_connection.cursor()

    @classmethod
    def tearDownClass(cls):
        cls.mock_resource.deactivate()
        # Removing test projection entries from projections db
        cls.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_genbank_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        driver = genbank.GenbankDriver('test')

        projection_configuration = PrototypeDeserializer('tests/test_genbank_config.yaml')
        self.genbank_projector = DBProjector('test_genbank_projection', driver,
                                             projection_configuration.prototype_tree,
                                             projection_configuration.root_projection_uri)

    def tearDown(self):
        # Clean up previous test entries
        self.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_genbank_projection' ")
        self.db_connection.commit()

    def test_create_projections(self):
        """
        Tests if Genbank projection manager creates projections
        """
        self.cursor.execute(" SELECT path FROM tree_table WHERE projection_name='test_genbank_projection' ")

        created_projections = [os.path.join(*r[0]) for r in self.cursor]

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(4, len(created_projections),
                         msg='Checking if Genbank projector created 4 projections, current number: {}'.format(len(created_projections)))

        # Check query projection creation
        self.assertIn('/test_query_folder',
                      created_projections,
                      msg='Checking creation of query projection')

        # Check gb file projection creation
        self.assertIn('/test_query_folder/sequence.gb',
                      created_projections,
                      msg='Checking creation of experiment projection.')

        # Check fasta projection creation
        self.assertIn('/test_query_folder/sequence.fasta',
                      created_projections,
                      msg='Checking creation of metadata projection.')


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
        cls.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_genbank_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        self.gbk_proj = subprocess.Popen([sys.executable,
                                          'genbank.py',
                                          '-p', 'test_genbank_projection',
                                          '-m', 'tests/mnt',
                                          '-d', 'tests/data',
                                          '-c', 'tests/test_genbank_config.yaml'],
                                         stdout=subprocess.DEVNULL)

        # Wait to initialize Projector
        time.sleep(0.5)

    def tearDown(self):
        # Shutting down genbank projection
        self.gbk_proj.terminate()

        # Unmounting projection dir
        subprocess.Popen(['fusermount', '-u', 'tests/mnt'])

        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_genbank_projection' ")
        self.db_connection.commit()

    def test_projections_contents(self):
        """
        Tests if file projections contents
        """
        # Loading test contents
        with open('tests/mock_resource/genbank_mock_data/mock_contents.txt') as f_f:
            test_contents = f_f.readlines()

        # Checking fasta file projection contents
        with open('tests/mnt/test_query_folder/sequence.fasta') as p_f:
            projected_fasta = p_f.readlines()

        self.assertEqual(projected_fasta, test_contents, msg='Check if fasta file contents loaded properly.')

        # Checking gb file projection contents
        with open('tests/mnt/test_query_folder/sequence.gb') as p_f:
            projected_gb = p_f.readlines()
        self.assertEqual(projected_gb, test_contents, msg='Check if gb file contents loaded properly.')
