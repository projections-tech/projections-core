import logging
import logging.config
from unittest import TestCase, skip
from tests.sra_mock import SRAMock
from projections import PrototypeDeserializer, Projection
import sra

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('sra_test')


class TestSRAProjectionManager(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mock_resource = SRAMock('http://eutils.ncbi.nlm.nih.gov', 'tests/mock_resource')

    def setUp(self):
        driver = sra.SRADriver('test')
        projection_configuration = PrototypeDeserializer('tests/test_sra_config.yaml')
        root_projection = Projection('/', projection_configuration.root_projection_uri)
        try:
            self.sra_projector = sra.SRAProjector(driver, root_projection, projection_configuration.prototype_tree)
        except:
            logger.debug('Encountered request to unmocked resource: %s', self.mock_resource.get_last_request_to_mock())

    def test_create_projections(self):
        """
        Tests if SRA projection manager creates projections
        """
        created_projections = self.sra_projector.projections

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(5, len(created_projections),
                         msg='Checking if SRA projector created 5 projections, current number: {}'.format(len(created_projections)))

        # Check query projection creation
        self.assertIn('/"Streptococcus"[Organism] OR Streptococcus[All Fields]',
                      created_projections,
                      msg='Checking creation of query projection')

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
