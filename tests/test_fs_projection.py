import getpass
import logging
import logging.config
import os
import subprocess
import sys
import time
from unittest import TestCase

import psycopg2

from db_projector import DBProjector
from drivers.fs_projection import FSDriver
from projections import PrototypeDeserializer

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'
CONFIG_PATH = 'tests/test_fs_proj_config.yaml'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('fs_proj_test')


class TestFSProjection(TestCase):
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
        cls.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_fs_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        driver = FSDriver()
        projection_configuration = PrototypeDeserializer(CONFIG_PATH)

        self.fs_projector = DBProjector('test_fs_projection', driver,
                                        projection_configuration.prototype_tree,
                                        projection_configuration.root_projection_uri)

    def tearDown(self):
        # Clean up previous test entries
        self.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_fs_projection' ")
        self.db_connection.commit()

    def test_create_projections(self):
        """
        Tests if filesystem projector creates projections
        """
        self.cursor.execute(" SELECT path FROM tree_table WHERE projection_name='test_fs_projection' ")

        created_projections = [os.path.join(*r[0]) for r in self.cursor]

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(37, len(created_projections),
                         msg='Checking if FS projector created 37 projections, current number: {0}'.format(len(created_projections)))

        # Check root dir projection creation
        self.assertIn('/test_dir',
                      created_projections,
                      msg='Checking creation of root projection')

        # Checking if FASTA files not present in config are not projected
        for i in range(1, 5):
            fasta_proj_path = '/test_dir/fasta_file_{0}.fasta'.format(i)
            self.assertNotIn(fasta_proj_path, created_projections,
                          msg='Checking if {0} not in projections'.format(fasta_proj_path))

        # Check BAM files projections creation
        for i in range(1, 6):
            bam_proj_path = '/test_dir/bam_file_{0}.bam'.format(i)
            self.assertIn(bam_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(bam_proj_path))

        for i in range(1, 6):
            # Check sample dir projection creation
            sample_dir_path = '/test_dir/sample_{0}'.format(i)
            self.assertIn(sample_dir_path, created_projections,
                          msg='Checking creation of {0} projection'.format(sample_dir_path))

            # Check sample VCF and BED projections creation
            vcf_proj_path = os.path.join(sample_dir_path, 'vcf_file.vcf')
            bed_proj_path = os.path.join(sample_dir_path, 'bed_file.bed')
            self.assertIn(vcf_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(vcf_proj_path))
            self.assertIn(bed_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(bed_proj_path))

            # Check inner rerun projection creation
            rerun_proj_path = os.path.join(sample_dir_path, 'rerun')
            self.assertIn(rerun_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(rerun_proj_path))

            # Check rerun contents BED and VCF files
            vcf_proj_path = os.path.join(rerun_proj_path, 'vcf_file.vcf')
            bed_proj_path = os.path.join(rerun_proj_path, 'bed_file.bed')
            self.assertIn(vcf_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(vcf_proj_path))
            self.assertIn(bed_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(bed_proj_path))


class TestFSProjectionContents(TestCase):
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
        cls.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_fs_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        self.fs_proj = subprocess.Popen([sys.executable,
                                         'fs_projection.py',
                                         '-p', 'test_fs_projection',
                                         '-m', 'tests/mnt',
                                         '-d', 'tests/data',
                                         '-c', 'tests/test_fs_proj_config.yaml'],
                                        stdout=subprocess.DEVNULL)

        # Wait to initialize Projector
        time.sleep(0.5)

    def tearDown(self):
        # Shutting down genbank projection
        self.fs_proj.terminate()

        # Unmounting projection dir
        subprocess.Popen(['fusermount', '-u', 'tests/mnt'])

        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_fs_projection' ")
        self.db_connection.commit()

    def test_projections_contents(self):
        # Unmount any previous tests
        subprocess.Popen(['fusermount', '-u', MOUNT_POINT])

        # Expected mock files contents
        expected_bam = 'Mock bam here!\n'
        expected_vcf = 'Mock vcf here!\n'
        expected_bed = 'Mock bed here!\n'
        for i in range(1, 6):
            # Paths to projections
            vcf_proj_path = 'tests/mnt/test_dir/sample_{}/vcf_file.vcf'.format(i)
            bam_proj_path = 'tests/mnt/test_dir/bam_file_{0}.bam'.format(i)
            bed_proj_path = 'tests/mnt/test_dir/sample_{}/bed_file.bed'.format(i)
            rerun_vcf_proj_path = 'tests/mnt/test_dir/sample_{}/rerun/vcf_file.vcf'.format(i)
            rerun_bed_proj_path = 'tests/mnt/test_dir/sample_{}/rerun/bed_file.bed'.format(i)

            # Reading projections contents on paths
            with open(bam_proj_path) as bam_f:
                test_bam_file = bam_f.read()

            with open(vcf_proj_path) as vcf_f:
                test_vcf_file = vcf_f.read()

            with open(bed_proj_path) as bed_f:
                test_bed_file = bed_f.read()

            with open(rerun_vcf_proj_path) as r_vcf_f:
                test_r_vcf_file = r_vcf_f.read()

            with open(rerun_bed_proj_path) as r_bed_f:
                test_r_bed_file = r_bed_f.read()

            # Checking projected files contents
            self.assertEqual(expected_bam, test_bam_file,
                             msg='Checking contents of projection {0}'.format(bam_proj_path))
            self.assertEqual(expected_vcf, test_vcf_file,
                             msg='Checking contents of projection {0}'.format(vcf_proj_path))
            self.assertEqual(expected_bed, test_bed_file,
                             msg='Checking contents of projection {0}'.format(bed_proj_path))
            self.assertEqual(expected_vcf, test_r_vcf_file,
                             msg='Checking contents of projection {0}'.format(rerun_vcf_proj_path))
            self.assertEqual(expected_bed, test_r_bed_file,
                             msg='Checking contents of projection {0}'.format(rerun_bed_proj_path))
