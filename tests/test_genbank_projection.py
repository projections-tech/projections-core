
import logging
import logging.config
from unittest import TestCase, skip
from tests.genbank_mock import GenbankMock
from projections import PrototypeDeserializer, Projection
import genbank
import subprocess
import httpretty
import time

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('genbank_test')


class TestGenbankProjector(TestCase):
    @classmethod
    def setUpClass(cls):
        httpretty.disable()
        cls.mock_resource = GenbankMock('http://eutils.ncbi.nlm.nih.gov', 'tests/mock_resource')

    def setUp(self):
        driver = genbank.GenbankDriver('test')
        projection_configuration = PrototypeDeserializer('tests/test_genbank_config.yaml')
        root_projection = Projection('/', projection_configuration.root_projection_uri)
        self.genbank_projector = genbank.GenbankProjector(driver, root_projection, projection_configuration.prototype_tree)

    def test_create_projections(self):
        """
        Tests if Genbank projection manager creates projections
        """
        created_projections = self.genbank_projector.projections

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(4, len(created_projections),
                         msg='Checking if Genbank projector created 4 projections, current number: {}'.format(len(created_projections)))

        # Check query projection creation
        self.assertIn('/GI:939732440',
                      created_projections,
                      msg='Checking creation of query projection')

        # Check gb file projection creation
        self.assertIn('/GI:939732440/sequence.gb',
                      created_projections,
                      msg='Checking creation of experiment projection.')

        # Check fasta projection creation
        self.assertIn('/GI:939732440/sequence.fasta',
                      created_projections,
                      msg='Checking creation of metadata projection.')

    def test_projections_contents(self):
        """
        Tests if file projections contents
        """
        # Starting Genbank projector process
        gbk_proj = subprocess.Popen(['./genbank.py',
                                     '-m', 'tests/mnt',
                                     '-d', 'tests/data',
                                     '-c', 'tests/test_genbank_config.yaml'],
                                    stdout=subprocess.DEVNULL)
        # Time to initialize Projector
        time.sleep(0.5)
        # Loading test contents
        with open('tests/mock_resource/genbank_mock_data/mock_contents.txt') as f_f:
            test_contents = f_f.readlines()
        # Checking fasta file projection contents
        with open('tests/mnt/GI:939732440/sequence.fasta') as p_f:
            projected_fasta = p_f.readlines()
        self.assertEqual(projected_fasta, test_contents, msg='Check if fasta file contents loaded properly.')
        # Checking gb file projection contents
        with open('tests/mnt/GI:939732440/sequence.gb') as p_f:
            projected_gb = p_f.readlines()
        self.assertEqual(projected_gb, test_contents, msg='Check if gb file contents loaded properly.')
        # Stopping projector
        gbk_proj.terminate()
        subprocess.Popen(['fusermount', '-u', 'tests/mnt'])
