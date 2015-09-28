__author__ = 'abragin'

import logging
import logging.config
from unittest import TestCase, skip

import iontorrent

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_iontorrent_projection')

# TODO: remove connection parameters to configuration
HOST = '10.5.20.17'
USER = 'ionadmin'
PASSWORD = '0ECu1lW'


class TestIonTorrentProjection(TestCase):

    def setUp(self):
        self.iontorrent = iontorrent.IonTorrentProjection(HOST, USER, PASSWORD)

    def test_create_projections(self):
        """
        Tests if ion torrent projection manager creates projections
        """
        self.iontorrent.create_projections()
        self.assertGreater(len(self.iontorrent.projections), 1)
        # TODO rewrite tests to use mock data