import os
import shutil
import time
import subprocess
import logging
import logging.config
from unittest import TestCase, skip

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('transparent_projections_test')


class TestTransparentProjection(TestCase):
    def setUp(self):
        # Clean data folder before tests
        for path in os.listdir(DATA_FOLDER):
            if path not in ['folder', 'file']:
                path = os.path.join(DATA_FOLDER, path)
                if os.path.isfile(path):
                    os.remove(path)
                else:
                    shutil.rmtree(path)

    def test_flat_iontorrent_projection(self):
        """
        Test transparent projection creation using Torrent Suite
        """
        subprocess.Popen(['fusermount', '-u', MOUNT_POINT])
        # Starting iontorrent projection
        ts_proj = subprocess.Popen(['./iontorrent.py',
                                    '-m', MOUNT_POINT,
                                    '-d', DATA_FOLDER,
                                    '-c', 'tests/test_ts_transparent_proj.yaml'],
                                   stdout=subprocess.DEVNULL)

        # Time to initialize projector properly
        time.sleep(0.5)

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

        # Terminating Projector process
        ts_proj.terminate()

        # Unmounting mnt dir
        subprocess.Popen(['fusermount', '-u', MOUNT_POINT])

    def test_flat_sra_projection(self):
        """
        Test transparent projection creation using SRA
        """
        subprocess.Popen(['fusermount', '-u', MOUNT_POINT])
        # Starting iontorrent projection
        sra_proj = subprocess.Popen(['./sra.py',
                                    '-m', MOUNT_POINT,
                                    '-d', DATA_FOLDER,
                                    '-c', 'tests/test_sra_transparent_proj.yaml'],
                                   stdout=subprocess.DEVNULL)

        # Time to initialize projector properly
        time.sleep(1)


        root_dir_contents = ['file', 'folder', 'SRX1058124_metadata.json', 'SRR2062160.sam']

        self.assertSetEqual(set(root_dir_contents), set(os.listdir('tests/mnt')), msg='Checking flat SRA projection.')

        # Terminating Projector process
        sra_proj.terminate()

        # Unmounting mnt dir
        subprocess.Popen(['fusermount', '-u', MOUNT_POINT])


