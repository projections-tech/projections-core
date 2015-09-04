__author__ = 'abragin'

import logging
import logging.config
from unittest import TestCase, skip

import iontorrent

# Import logging configuration from the file provided
logging.config.fileConfig('../logging.cfg')
logger = logging.getLogger('test_filesystem')

HOST = '10.5.20.17'
USER = 'ionadmin'
PASSWORD = '0ECu1lW'


class TestIonTorrentProjection(TestCase):

    def setUp(self):
        self.iontorrent = iontorrent.IonTorrentProjection(HOST, USER, PASSWORD)

    def test_create_projections(self):
        self.iontorrent.create_projections()