__author__ = 'abragin'

import os
import httpretty
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

    def test_create_projections(self):
        """
        Tests if ion torrent projection manager creates projections
        """
        self.iontorrent.create_projections()
        self.assertGreater(len(self.iontorrent.projections), 1)

    def test_bam_projections_creation(self):
        self.iontorrent.create_projections()
        bam_files_list = [p for p in self.iontorrent.projections if os.path.splitext(p)[1] == '.bam']
        logger.debug('Bam files projections: %s', bam_files_list)
        logger.debug('Mock Torrent Suite requests: %s', self.mock_resource.get_last_request_to_mock())
        self.assertEqual(len(bam_files_list), 5)

    def test_vcf_projections_creation(self):
        self.iontorrent.create_projections()
        vcf_files_list = [p for p in self.iontorrent.projections if os.path.splitext(p)[1] == '.vcf']
        self.assertGreater(len(vcf_files_list), 1)
