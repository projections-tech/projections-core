__author__ = 'abragin'

import logging
import logging.config
import os
import time
from subprocess import Popen
from unittest import TestCase, skip

import filesystem

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_filesystem')

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

class TestFilesystem(TestCase):

    def setUp(self):
        logger.info('Mounting Projections filesystem for testing')
        Popen(['./filesystem.py', MOUNT_POINT, DATA_FOLDER])
        # Let filesystem time to mount
        time.sleep(0.5)

    def tearDown(self):
        logger.info('Unmounting testing Projections filesystem')
        Popen(['fusermount', '-u', MOUNT_POINT])

    def test_mounting(self):
        mount_content = os.listdir(MOUNT_POINT)
        logger.debug('Mounted directory content: %s', mount_content)

        # Check that file data folder content is visible via mount point
        self.assertTrue('folder' in mount_content)
        self.assertTrue('file' in mount_content)
        dir_content = os.listdir(os.path.join(MOUNT_POINT, 'folder'))
        self.assertTrue('file' in dir_content)

    def test_hello_world_projection(self):
        logger.debug('Testing hello world projection')
        mount_content = os.listdir(MOUNT_POINT)
        self.assertTrue('projection' in mount_content, 'Projection resource is present in mount root')
        dir_content = os.listdir(os.path.join(MOUNT_POINT, 'folder'))
        self.assertFalse('projection' in dir_content, 'Projection is not present below first level of mount folder')

        projection_path = os.path.join(MOUNT_POINT, 'projection')

        self.assertEqual(1, os.stat(projection_path).st_size, 'Check projection size before reading')

        with open(projection_path) as f:
            content = f.read()

        self.assertEqual('Hello World!\n', content, 'Check projection content')
        self.assertEqual(len('Hello World!\n'), os.stat(projection_path).st_size, 'Check projection size after reading')

        self.assertTrue('projection' in os.listdir(DATA_FOLDER), 'Check that file now in data folder')
        self.assertEquals(1, len([f for f in os.listdir(DATA_FOLDER) if f == 'projection']), 'Check that projection is present only once')

        # Remove projection resource
        os.remove(projection_path)

        self.assertFalse('projection' in os.listdir(DATA_FOLDER), 'Check that file is not more in data folder')
        self.assertEqual(1, os.stat(projection_path).st_size, 'Check projection size after deletion')
