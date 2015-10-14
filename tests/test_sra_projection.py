import os
import logging
import logging.config
from unittest import TestCase, skip
from .mock import SRAMock


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
        os.environ["http_proxy"] = "http://mocksra.com"
        cls.mock_resource = SRAMock('mocksra.com', 'tests/mock_resource')

    def setUp(self):
        self.mock_url = self.mock_resource.mock_url
        driver = sra.SRADriver('vsvekolkin@parseq.pro')
        try:
            self.sra_projector = sra.SRAProjector(driver)
        except:
            logger.debug('Last request to mock: %s', self.mock_resource.get_last_request_to_mock())

    def test_create_projections(self):
        """
        Tests if SRA projection manager creates projections
        """
        self.assertGreater(len(self.sra_projector.projections), 1)

    def test_sam_projections_creation(self):
        sam_files_list = [p for p in self.sra_projector.projections if os.path.splitext(p)[1] == '.sam']
        logger.debug('Bam files projections: %s', sam_files_list)
        self.assertEqual(len(sam_files_list), 5)
