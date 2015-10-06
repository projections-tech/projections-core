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
        expected_experiment_path = '/Sequoia SN1-14-Newborn 5 samples from CV'
        projection_paths = [key for key in self.iontorrent.projections]
        paths_experiment_dir_state = all([re.search(expected_experiment_path, path) for path in projection_paths])
        self.assertTrue(paths_experiment_dir_state)

    def test_run_projection_creattion(self):
        expected_run_name = 'Run_11_hg19_v3_008'
        projection_paths = [key for key in self.iontorrent.projections]
        paths_run_dir_state = all([re.search(expected_run_name, path) for path in projection_paths])
        self.assertTrue(paths_run_dir_state, paths_run_dir_state)

    def test_bam_projections_creation(self):
        bam_files_list = sorted([os.path.basename(p) for p in self.iontorrent.projections if os.path.splitext(p)[1] == '.bam'])
        expected_bam_files = ['sample_{}.bam'.format(i+1) for i in range(5)]
        self.assertListEqual(bam_files_list, expected_bam_files)

    def test_vcf_projections_creation(self):
        vcf_files_list = [p for p in self.iontorrent.projections if os.path.splitext(p)[1] == '.vcf']
        self.assertGreater(len(vcf_files_list), 1)
