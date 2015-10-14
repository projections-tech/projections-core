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
        created_projections = self.sra_projector.projections

        logger.debug(created_projections.keys())

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(5, len(created_projections),
                         msg='Checking if SRA projector created 5 projections, current number: {}'.format(len(created_projections)))

        # Check query projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]', created_projections,
                      msg='Checking creation of query projection')

        # Check experiment projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]/SRX1058124', created_projections,
                      msg='Checking creation of experiment projection.')

        # Check metadata projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]/SRX1058124/metadata.json', created_projections,
                      msg='Checking creation of metadata projection.')

        # Check sam file projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]/SRX1058124/SRR2062160.sam', created_projections,
                      msg='Checking creation of sam file projection.')
