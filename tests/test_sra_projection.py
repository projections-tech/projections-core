import logging
import logging.config
from unittest import TestCase, skip
from tests.sra_mock import SRAMock
from projections import PrototypeDeserializer, Projector
import sra

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('sra_test')


class SRAProjectionManager(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mock_resource = SRAMock('http://eutils.ncbi.nlm.nih.gov', 'tests/mock_resource')

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
        self.assertEqual(4, len(created_projections),
                         msg='Checking if SRA projector created 4 projections, current number: {}'.format(len(created_projections)))

        # Check experiment projection creation
        self.assertIn('/SRX1058124',
                      created_projections,
                      msg='Checking creation of experiment projection.')

        # Check metadata projection creation
        self.assertIn('/SRX1058124/metadata.json',
                      created_projections,
                      msg='Checking creation of metadata projection.')

        # Check sam file projection creation
        self.assertIn('/SRX1058124/SRR2062160.sam',
                      created_projections,
                      msg='Checking creation of sam file projection.')
