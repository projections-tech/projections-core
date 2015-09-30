__author__ = 'abragin'

import os
import time
import subprocess
import logging
import logging.config
from fuse import FUSE
from unittest import TestCase, skip
from .resource_pretender import Torrent_Suite_Mock



import iontorrent
from filesystem import ProjectionFilesystem

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_iontorrent_projection')

# TODO: remove connection parameters to configuration
USER = 'ionadmin'
PASSWORD = '0ECu1lW'


class TestIonTorrentProjection(TestCase):
    def test_create_projections(self):
        """
        Tests if ion torrent projection manager creates projections
        """
        self.mock_resource = Torrent_Suite_Mock('tests/mock_resource')
        self.mock_url = self.mock_resource.mock_url
        self.iontorrent = iontorrent.IonTorrentProjection(self.mock_url, USER, PASSWORD)
        self.iontorrent.create_projections()
        self.assertGreater(len(self.iontorrent.projections), 1)

    def test_bam_files_creation(self):
        self.mock_resource = Torrent_Suite_Mock('tests/mock_resource')
        self.mock_url = self.mock_resource.mock_url
        self.iontorrent = iontorrent.IonTorrentProjection(self.mock_url, USER, PASSWORD)
        logger.debug('Current mock server: %s', self.mock_url)
        self.iontorrent.create_projections()
        projection_filesystem = ProjectionFilesystem(MOUNT_POINT, DATA_FOLDER)
        projection_filesystem.projection_manager = self.iontorrent
        fuse = FUSE(projection_filesystem, MOUNT_POINT, foreground=True, nonempty=True)
        logger.debug(os.listdir(MOUNT_POINT))
