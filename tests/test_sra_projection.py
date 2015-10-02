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
logger = logging.getLogger('test_sra_projection')

# TODO: remove connection parameters to configuration
USER = 'ionadmin'
PASSWORD = '0ECu1lW'

class SRAProjectionManager(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mock_resource = SRAMock('mocksra.com','tests/mock_resource')

    def setUp(self):
        self.mock_url = self.mock_resource.mock_url
        self.sra = sra.SRAProjectionManager('vsvekolkin@parseq.pro', 'Streptococcus', 1)

    def test_create_projections(self):
        """
        Tests if SRA projection manager creates projections
        """
        self.sra.create_projections()
        self.assertGreater(len(self.sra.projections), 1)

    def test_sam_projections_creation(self):
        self.sra.create_projections()
        sam_files_list = [p for p in self.sra.projections if os.path.splitext(p)[1] == '.sam']
        logger.debug('Bam files projections: %s', sam_files_list)
        self.assertEqual(len(sam_files_list), 5)
