__author__ = 'abragin'

import logging
import logging.config
from unittest import TestCase, skip
from .resource_pretender import Torrent_Suite_Mock


import iontorrent

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_iontorrent_projection')

# TODO: remove connection parameters to configuration
USER = 'ionadmin'
PASSWORD = '0ECu1lW'


class TestIonTorrentProjection(TestCase):
    def setUp(self):
        self.mock_resource = Torrent_Suite_Mock('tests/mock_resource')
        self.mock_url = self.mock_resource.mock_url
        self.iontorrent = iontorrent.IonTorrentProjection(self.mock_url, USER, PASSWORD)


    def test_create_projections(self):
        """
        Tests if ion torrent projection manager creates projections
        """
        self.iontorrent.create_projections()
        logger.debug('Torrent Suite projections: %s', self.iontorrent.projections.items())
        self.assertGreater(len(self.iontorrent.projections), 1)