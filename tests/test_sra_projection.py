import getpass
import logging
import logging.config
import os
from unittest import TestCase

import psycopg2

import drivers.sra_driver as sra
from db_projector import DBProjector
from projections import PrototypeDeserializer
from tests.mock import MockResource

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('sra_test')


class TestSRAProjector(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.mock_resource = MockResource('tests/sra_mock.json')

        # Initializing database connection which will be used during tests
        cls.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))
        # Creating cursor, which will be used to interact with database
        cls.cursor = cls.db_connection.cursor()

    @classmethod
    def tearDownClass(cls):
        cls.mock_resource.deactivate()

        # Removing test projection entries from projections db
        cls.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_sra_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        driver = sra.SRADriver('test')
        projection_configuration = PrototypeDeserializer('tests/test_sra_config.yaml')

        self.sra_projector = DBProjector('test_sra_projection', driver,
                                         projection_configuration.prototype_tree,
                                         projection_configuration.root_projection_uri)

    def tearDown(self):
        # Clean up previous test entries
        self.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_sra_projection' ")
        self.db_connection.commit()

    def test_create_projections(self):
        """
        Tests if SRA projection manager creates projections
        """
        self.cursor.execute(" SELECT path FROM tree_table WHERE projection_name='test_sra_projection' ")

        created_projections = [os.path.join(*r[0]) for r in self.cursor]

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(5, len(created_projections),
                         msg='SRA projector created {} projections.'.format(len(created_projections)))

        # Check search query projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]',
                      created_projections,
                      msg='Checking creation of search query projection.')

        # Check experiment projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]/SRX1058124',
                      created_projections,
                      msg='Checking creation of experiment projection.')

        # Check metadata projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]/SRX1058124/metadata.json',
                      created_projections,
                      msg='Checking creation of metadata projection.')

        # Check sam file projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]/SRX1058124/SRR2062160.sam',
                      created_projections,
                      msg='Checking creation of sam file projection.')
