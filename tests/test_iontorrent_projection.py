__author__ = 'abragin'

import logging
import logging.config
import os
import time
from unittest import TestCase, skip
from subprocess import Popen

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

    def tearDown(self):
        Popen(['fusermount', '-u', MOUNT_POINT])

    def test_create_projections(self):
        self.iontorrent.create_projections()
        self.assertGreater(len(self.iontorrent.projections), 1)