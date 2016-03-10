import getpass
import logging
import logging.config
import os
import shutil
import subprocess
import sys
import time
from unittest import TestCase

import psycopg2

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('transparent_projections_test')


class TestTransparentIontorrentProjection(TestCase):
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
        cls.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_iontorrent_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        # Clean data folder before tests
        for path in os.listdir(DATA_FOLDER):
            if path not in ['folder', 'file']:
                path = os.path.join(DATA_FOLDER, path)
                if os.path.isfile(path):
                    os.remove(path)
                else:
                    shutil.rmtree(path)

        self.iontorrent_projection = subprocess.Popen([sys.executable,
                                                       'iontorrent.py',
                                                       '-p', 'test_iontorrent_projection',
                                                       '-m', 'tests/mnt',
                                                       '-d', 'tests/data',
                                                       '-c', 'tests/test_ts_transparent_proj.yaml'],
                                                      stdout=subprocess.DEVNULL)

        # Wait to initialize Projector
        time.sleep(1)

    def tearDown(self):
        # Shutting down genbank projection
        self.iontorrent_projection.terminate()
        logger.debug('Unmounting mount dir!')
        # Unmounting projection dir
        subprocess.Popen(['fusermount', '-u', 'tests/mnt'])

        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_iontorrent_projection' ")
        self.db_connection.commit()

    def test_flat_iontorrent_projection(self):
        """
        Test transparent projection creation using Torrent Suite
        """
        root_dir_contents = ['IAD66589_181_SNP-HID-p2-L_Target_regions.bed', 'file', 'folder',
                             'test_experiment_1_metadata.json', 'test_experiment_1_plannedexperiment.json',
                             'test_experiment_1_sample_1_TSVC_variants.vcf',
                             'test_experiment_1_sample_2_TSVC_variants.vcf', 'test_experiment_2_metadata.json',
                             'test_experiment_2_plannedexperiment.json', 'test_experiment_2_sample_3_TSVC_variants.vcf',
                             'test_experiment_2_sample_4_TSVC_variants.vcf', 'test_run_1_sample_1.bam',
                             'test_run_1_sample_1_metadata.json', 'test_run_1_sample_2.bam',
                             'test_run_1_sample_2_metadata.json', 'test_run_1_variant_caller_settings.json',
                             'test_run_2_sample_1.bam', 'test_run_2_sample_1_metadata.json', 'test_run_2_sample_2.bam',
                             'test_run_2_sample_2_metadata.json', 'test_run_2_variant_caller_settings.json',
                             'test_run_3_sample_3.bam', 'test_run_3_sample_3_metadata.json', 'test_run_3_sample_4.bam',
                             'test_run_3_sample_4_metadata.json', 'test_run_3_variant_caller_settings.json',
                             'test_run_4_sample_3.bam', 'test_run_4_sample_3_metadata.json', 'test_run_4_sample_4.bam',
                             'test_run_4_sample_4_metadata.json', 'test_run_4_variant_caller_settings.json']

        self.assertSetEqual(set(root_dir_contents), set(os.listdir('tests/mnt')),
                            msg='Checking flat iontorrent projection.')


class TestTransparentSRAProjection(TestCase):
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
        cls.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_sra_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        # Clean data folder before tests
        for path in os.listdir(DATA_FOLDER):
            if path not in ['folder', 'file']:
                path = os.path.join(DATA_FOLDER, path)
                if os.path.isfile(path):
                    os.remove(path)
                else:
                    shutil.rmtree(path)

        self.sra_projection = subprocess.Popen([sys.executable,
                                                'sra.py',
                                                '-p', 'test_sra_projection',
                                                '-m', 'tests/mnt',
                                                '-d', 'tests/data',
                                                '-c', 'tests/test_sra_transparent_proj.yaml'],
                                               stdout=subprocess.DEVNULL)

        # Wait to initialize Projector
        time.sleep(1)

    def tearDown(self):
        # Shutting down genbank projection
        self.sra_projection.terminate()

        logger.debug('Unmounting mount dir!')
        # Unmounting projection dir
        subprocess.Popen(['fusermount', '-u', 'tests/mnt'])

        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_sra_projection' ")
        self.db_connection.commit()

    def flat_sra_projection(self):
        """
        Test transparent projection creation using SRA
        """

        root_dir_contents = ['file', 'folder', 'metadata.json', 'SRR2062160.sam']

        self.assertSetEqual(set(root_dir_contents), set(os.listdir('tests/mnt')), msg='Checking flat SRA projection.')