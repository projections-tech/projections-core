import os
import logging
import logging.config
from unittest import TestCase, skip
from .mock import SRAMock
import httpretty


import sra

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('sra_test')

# TODO: remove connection parameters to configuration
USER = 'ionadmin'
PASSWORD = '0ECu1lW'

class SRAProjectionManager(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mock_resource = SRAMock('http://eutils.ncbi.nlm.nih.gov', 'tests/mock_resource')

    def setUp(self):
        self.mock_url = self.mock_resource.mock_url
        driver = sra.SRADriver('test')
        try:
            self.sra_projector = sra.SRAProjector(driver)
        except:
            logger.debug(self.mock_resource.get_last_request_to_mock())

    def test_create_projections(self):
        """
        Tests if SRA projection manager creates projections
        """
        self.assertGreater(len(self.sra_projector.projections), 1)
