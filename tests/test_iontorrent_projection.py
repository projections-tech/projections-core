__author__ = 'abragin'

import os
import re
import logging
import logging.config
from unittest import TestCase, skip
from .mock import TorrentSuiteMock



import iontorrent

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_iontorrent_projection')

USER = 'user'
PASSWORD = 'password'

class TestIonTorrentProjection(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mock_resource = TorrentSuiteMock('mockiontorrent.com','tests/mock_resource')
        cls.expected_experiment_path = '/Sequoia SN1-14-Newborn 5 samples from CV'
        cls.expected_run_name = 'Run_11_hg19_v3_008'

    def setUp(self):
        self.mock_url = self.mock_resource.mock_url
        self.iontorrent = iontorrent.IonTorrentProjection(self.mock_url, USER, PASSWORD)
        self.iontorrent.create_projections()

    def test_create_projections(self):
        """
        Tests if ion torrent projection manager creates projections
        """

        self.assertGreater(len(self.iontorrent.projections), 1)

    def test_experiment_projection_creation(self):
        projection_paths = [key for key in self.iontorrent.projections]
        paths_experiment_dir_state = all([re.search(self.expected_experiment_path, path) for path in projection_paths])
        self.assertTrue(paths_experiment_dir_state)

    def test_run_projection_creation(self):
        # All projections paths except '/Sequoia SN1-14-Newborn 5 samples from CV' and its metadata files
        projection_paths = [key for key in self.iontorrent.projections if not re.search(self.expected_experiment_path+'$', key)
                            and re.search(self.expected_experiment_path + '/plannedexperiment.json$', key)
                            and re.search(self.expected_experiment_path + '/metadata.json$', key)]
        paths_run_dir_state = all([re.search(self.expected_run_name, path) for path in projection_paths])
        self.assertTrue(paths_run_dir_state)

    def test_bam_projections_creation(self):
        bam_files_list = sorted([os.path.basename(p) for p in self.iontorrent.projections if os.path.splitext(p)[1] == '.bam'])
        expected_bam_files = ['sample_{}.bam'.format(i) for i in range(1,6)]
        self.assertListEqual(bam_files_list, expected_bam_files)

    def test_vcf_projections_creation(self):
        vcf_files_list = sorted([p for p in self.iontorrent.projections if os.path.splitext(p)[1] == '.vcf'])
        for sample_id in range(1,6):
            for variant_caller_dir in ['', '.49', '.50']:
                for variant_file_name in ['TSVC_variants.vcf', 'all.merged.vcf', 'indel_assembly.vcf',
                                                              'indel_variants.vcf', 'small_variants.left.vcf',
                                                              'small_variants.vcf', 'small_variants_filtered.vcf',
                                                              'small_variants.sorted.vcf', 'SNP_variants.vcf']:
                    vcf_file_path = os.path.join(self.expected_experiment_path,
                                                 self.expected_run_name,
                                                 'sample_{}'.format(sample_id),
                                                 'variantCaller_out{}'.format(variant_caller_dir),
                                                 variant_file_name)
                    self.assertTrue(vcf_file_path in vcf_files_list, msg='VCF path: {0}'.format(vcf_file_path))