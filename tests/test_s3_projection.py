import logging
import logging.config
from unittest import TestCase, skip
from projections import PrototypeDeserializer, Projection
from aws_s3 import S3Driver, S3Projector

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('s3_test')

KEY_ID = 'AKIAIONUXTO6TR3UU3TQ'
ACCESS_KEY = 'UgYV9YRRoFX64nmxoL+4ry3QLBD0rPdoQRVTCB5w'
REGION_NAME = 'us-west-2'


class TestS3Driver(TestCase):

    def setUp(self):
        self.driver = S3Driver(KEY_ID, ACCESS_KEY, REGION_NAME, 'parseq')

    def test_bucket_contents(self):
        """
        Tests if driver gets proper bucket on init
        """
        exp_contents = ['projects/', 'projects/ensembl.txt']
        for obj in exp_contents:
            self.assertIn(obj, [o.key for o in self.driver.bucket.objects.all()], msg='Checking objects existence.')

    def test_get_uri_contents_as_stream(self):
        """
        Tests driver get_uri_contents_stream method
        """
        with open('tests/test_ensembl.txt', 'rb') as t_f:
            test_ensembl = t_f.read()
        for key, exp_content in {'projects/': b'', 'projects/ensembl.txt': test_ensembl}.items():
            self.assertEqual(exp_content, self.driver.get_uri_contents_as_stream(key),
                             msg='Checking content of object with URI: {0}'.format(key))

    def test_get_uri_contents_as_dict(self):
        """
        Tests driver get_uri_contents_as_dict method
        """
        contents_meta = [{'content_encoding': None, 'size': 0, 'name': 'projects/',
                          'content_type': 'binary/octet-stream', 'metadata': {}},
                         {'content_encoding': None, 'size': 1387, 'name': 'projects/ensembl.txt',
                          'content_type': 'text/plain', 'metadata': {'madefor': 'testing', 'quality': 'good'}}]
        for meta in contents_meta:
            self.assertDictEqual(meta, self.driver.get_uri_contents_as_dict(meta['name']),
                                 msg='Checking meta of object with URI: {0}'.format(meta['name']))


class TestS3Projector(TestCase):
    def setUp(self):
        projection_configuration = PrototypeDeserializer('tests/test_s3.yaml')

        root_projection = Projection('/', projection_configuration.root_projection_uri)

        s3_driver = S3Driver(KEY_ID, ACCESS_KEY, REGION_NAME, projection_configuration.root_projection_uri)
        self.s3_projector = S3Projector(s3_driver, root_projection, projection_configuration.prototype_tree)

    def test_create_projections(self):
        """
        Tests if S3 projector creates projections
        """
        created_projections = self.s3_projector.projections

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(6, len(created_projections),
                         msg='Checking if S3 projector created 6 projections, current number: {}'.format(len(created_projections)))

        # Check bucket projection creation
        self.assertIn('/parseq',
                      created_projections,
                      msg='Checking creation of bucket projection')

        # Check data file projection creation
        for d_f in ['/parseq/projects_ensembl.txt', '/parseq/projects_ensembl_1.txt',
                    '/parseq/ensembl_3.txt', '/parseq/ensembl_4.txt']:
            self.assertIn(d_f,
                          created_projections,
                          msg='Checking creation of {0} projection.'.format(d_f))

        # Check if files filtered by tag properly
        for d_f in ['/parseq/projects_ensembl_2.txt', '/parseq/projects_ensembl_3.txt',
                    '/parseq/projects_ensembl_4.txt', '/parseq/ensembl_1.txt', '/parseq/ensembl_2.txt']:
            self.assertNotIn(d_f, created_projections,
                             msg='Checking if projection {0} filtered.'.format(d_f))