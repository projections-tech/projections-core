import getpass
import logging
import logging.config
import os
from unittest import TestCase

import psycopg2

import fs_projection
from db_projector import DBProjector
from projections import PrototypeDeserializer

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'
CONFIG_PATH = 'tests/test_miseq_config.yaml'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('miseq_test')


class TestMiseqProjection(TestCase):
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
        cls.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_miseq_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        driver = fs_projection.FSDriver()
        projection_configuration = PrototypeDeserializer(CONFIG_PATH)

        self.miseq_projector = DBProjector('test_miseq_projection', driver,
                                           projection_configuration.prototype_tree,
                                           projection_configuration.root_projection_uri)

    def tearDown(self):
        # Clean up previous test entries
        self.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_miseq_projection' ")
        self.db_connection.commit()

    def test_create_projections(self):
        """
        Tests if MiSeq projector creates projections
        """
        self.cursor.execute(" SELECT path FROM tree_table WHERE projection_name='test_miseq_projection' ")

        created_projections = [os.path.join(*r[0]) for r in self.cursor]

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(38, len(created_projections),
                         msg='MiSeq projector created {0} projections.'.format(len(created_projections)))

        # Check root projection creation
        self.assertIn('/',
                      created_projections,
                      msg='Checking creation of root projection.')

        # Check data dir projection creation
        self.assertIn('/data',
                      created_projections,
                      msg='Checking creation of data dir projection.')

        base_dir = '/data'
        for run_dir_name in ['151028_M70294_0001_000000000-AGLPP',
                             '151102_M70294_0002_000000000-AGLP0',
                             '151116_M70294_0003_000000000-AGLWG']:
            run_dir_path = os.path.join(base_dir, run_dir_name)
            self.assertIn(run_dir_path, created_projections,
                          msg='Checking creation of run dir projection.')

            sample_sheet_path = os.path.join(run_dir_path, 'SampleSheetUsed.csv')
            self.assertIn(sample_sheet_path, created_projections,
                          msg='Checking creation of Sample Sheet projection.')

            # 'Data', 'Intensities' and 'Alignment' dirs are transparent projections, checking if they are not projected
            data_dir_path = os.path.join(run_dir_path, 'Data')
            self.assertNotIn(data_dir_path, created_projections,
                             msg='Checking creation of Data dir projection.')

            intensities_dir_path = os.path.join(data_dir_path, 'Intensities')
            self.assertNotIn(intensities_dir_path, created_projections,
                             msg='Checking creation of Intensities dir projection.')

            alignment_dir_path = os.path.join(intensities_dir_path, 'Alignment')
            self.assertNotIn(alignment_dir_path, created_projections,
                             msg='Checking creation of Alignment dir projection.')

            # Checking bam and bai projection creation
            sample_name = run_dir_name.split('-')[-1]
            for i in range(1, 6):
                bam_file_path = os.path.join(run_dir_path, '{0}-{1}-{1}_S{1}.bam'.format(sample_name, i))
                bai_file_path = os.path.join(run_dir_path, '{0}-{1}-{1}_S{1}.bam.bai'.format(sample_name, i))

                self.assertIn(bam_file_path, created_projections,
                              msg='Checking creation of BAM file projection.')
                self.assertIn(bai_file_path, created_projections,
                              msg='Checking creation of BAI filke projection.')

    def test_projections_contents(self):
        """
        Test if ProjectionTree open_resource method returns correct file contents
        """
        base_dir = '/data'
        for run_dir_name in ['151028_M70294_0001_000000000-AGLPP',
                     '151102_M70294_0002_000000000-AGLP0',
                     '151116_M70294_0003_000000000-AGLWG']:
            proj_run_dir_path = os.path.join(base_dir, run_dir_name)

            # Checking SampleSheetUsed.csv contents
            sample_sheet_path = os.path.join('tests', 'mock_resource',
                                             'miseq_mock_data', 'data', run_dir_name, 'Data',
                                             'Intensities', 'BaseCalls', 'Alignment', 'SampleSheetUsed.csv')
            with open(sample_sheet_path, 'rb') as ssc:
                sample_sheet_contents = ssc.read()

            proj_sample_sheet_path = os.path.join(proj_run_dir_path, 'SampleSheetUsed.csv')
            proj_sample_sheet_contents = self.miseq_projector.open_resource(proj_sample_sheet_path)[1].getvalue()

            self.assertEqual(sample_sheet_contents,
                             proj_sample_sheet_contents,
                             msg='Checking sample sheet contents correctness.')

            sample_name = run_dir_name.split('-')[-1]
            for i in range(1, 6):
                bam_file_path = os.path.join(proj_run_dir_path, '{0}-{1}-{1}_S{1}.bam'.format(sample_name, i))
                bai_file_path = os.path.join(proj_run_dir_path, '{0}-{1}-{1}_S{1}.bam.bai'.format(sample_name, i))

                # Checking if ProjectionTree gets correct BAM and BAI files contents on path
                proj_bam_contents = self.miseq_projector.open_resource(bam_file_path)[1].getvalue()
                self.assertEqual(proj_bam_contents,
                                 b'Mock bam here!\n',
                                 msg='Checking BAM file contents.')

                proj_bai_contents = self.miseq_projector.open_resource(bai_file_path)[1].getvalue()
                self.assertEqual(proj_bai_contents,
                                 b'Mock bai here!\n',
                                 msg='Checking BAM file contents.')







