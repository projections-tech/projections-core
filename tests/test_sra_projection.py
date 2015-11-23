import logging
import logging.config
from unittest import TestCase, skip
from tests.mock import MockResource
from projections import PrototypeDeserializer, Projector

import sra

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('sra_test')


class TestSRAProjector(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.mock_resource = MockResource('tests/sra_mock.json')

    @classmethod
    def tearDownClass(cls):
        cls.mock_resource.deactivate()

    def setUp(self):
        driver = sra.SRADriver('test')
        projection_configuration = PrototypeDeserializer('tests/test_sra_config.yaml')
        self.sra_projector = Projector(driver, projection_configuration.root_projection_uri,
                                       projection_configuration.prototype_tree)

    def test_create_projections(self):
        """
        Tests if SRA projection manager creates projections
        """
        created_projections = [n.get_path() for n in self.sra_projector.projection_tree.get_tree_nodes()]

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(5, len(created_projections),
                         msg='SRA projector created {} projections.'.format(len(created_projections)))

        # Check search query projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]',
                      created_projections,
                      msg='Checking creation of search query projection.')

        # Check experiment projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]/SRX1058124',
                      created_projections,
                      msg='Checking creation of experiment projection.')

        # Check metadata projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]/SRX1058124/metadata.json',
                      created_projections,
                      msg='Checking creation of metadata projection.')

        # Check sam file projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]/SRX1058124/SRR2062160.sam',
                      created_projections,
                      msg='Checking creation of sam file projection.')
