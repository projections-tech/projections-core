__author__ = 'abragin'

import logging
import logging.config
import os
from subprocess import Popen
from unittest import TestCase

# Import logging configuration from the file provided
logging.config.fileConfig('../logging.cfg')
logger = logging.getLogger('test_filesystem')

MOUNT_POINT = 'mnt'
DATA_FOLDER = 'data'

class TestFilesystem(TestCase):

    def setUp(self):
        logger.info('Mounting Projections filesystem for testing')
        self.mount_process = Popen(['../filesystem.py', MOUNT_POINT])

    def tearDown(self):
        logger.info('Unmounting testing Projections filesystem')
        Popen(['fusermount', '-u', MOUNT_POINT])
        self.mount_process.terminate()

    def test_mounting(self):
        mount_content = os.listdir(MOUNT_POINT)
        logger.debug('Mounted directory content: %s', mount_content)

        # Check that file data folder content is visible via mount point
        self.assertTrue('folder' in mount_content)
        self.assertTrue('file' in mount_content)
        dir_content = os.listdir(os.path.join(MOUNT_POINT, 'folder'))
        self.assertTrue('file' in dir_content)

    def test_hello_world_projection(self):
        mount_content = os.listdir(MOUNT_POINT)
        self.assertTrue('projection' in mount_content, 'Projection resource is present in mount root')
        dir_content = os.listdir(os.path.join(MOUNT_POINT, 'folder'))
        self.assertFalse('projection' in dir_content, 'Projection is not present below first level of mount folder')