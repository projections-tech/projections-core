import logging
import logging.config
from unittest import TestCase, skip
from projections import PrototypeDeserializer, Projector
from aws_s3 import S3Driver
from moto import mock_s3
import boto3
# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('s3_test')


class TestS3Driver(TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set`s up mock S3 resource
        """
        # Init S3 mock
        cls.mock = mock_s3()
        cls.mock.start()

        # Instantiate S3Driver which will be tested
        cls.driver = S3Driver('test_id', 'test_key', 'us-west-2', 'parseq')

        # Add contents to mock S3 resource using boto3
        s3 = boto3.resource('s3')
        s3.create_bucket(Bucket='parseq', CreateBucketConfiguration={'LocationConstraint': 'us-west-2'})
        s3.Object('parseq', 'projects/').put(Body=b'')
        s3.Object('parseq', 'projects/ensembl.txt').put(Body=b'Test ensembl here!',
                                                        Metadata={'madefor': 'testing', 'quality': 'good'})

    @classmethod
    def tearDownClass(cls):
        # Stopping mock resource
        cls.mock.stop()

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
        test_ensembl = b'Test ensembl here!'
        for key, exp_content in {'projects/': b'', 'projects/ensembl.txt': test_ensembl}.items():
            self.assertEqual(exp_content, self.driver.get_uri_contents_as_bytes(key),
                             msg='Checking content of object with URI: {0}'.format(key))

    def test_get_uri_contents_as_dict(self):
        """
        Tests driver get_uri_contents_as_dict method
        """
        contents_meta = [{'content_encoding': None, 'size': 0, 'name': 'projects/',
                          'content_type': 'text/plain; charset=utf-8', 'metadata': {}, 'resource_uri': 'projects/'},
                         {'content_encoding': None, 'size': 18, 'name': 'projects/ensembl.txt',
                          'content_type': 'text/plain; charset=utf-8',
                          'resource_uri': 'projects/ensembl.txt', 'metadata': {'madefor': 'testing', 'quality': 'good'}}]
        for meta in contents_meta:
            self.assertDictEqual(meta, self.driver.get_uri_contents_as_dict(meta['name']),
                                 msg='Checking meta of object with URI: {0}'.format(meta['name']))


class TestS3Projector(TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set`s up mock S3 resource
        """
        # Init S3 mock
        cls.mock = mock_s3()
        cls.mock.start()
        # Add contents to S3 resource using boto3
        s3 = boto3.resource('s3')
        s3.create_bucket(Bucket='parseq', CreateBucketConfiguration={'LocationConstraint': 'us-west-2'})
        s3.Object('parseq', 'projects/').put(Body=b'')

        # Setting 'quality' metadata field of files projections and adding files to mock resource
        for i in range(1, 5):
            if i <= 2:
                quality = 'bad'
            else:
                quality = 'good'
            s3.Object('parseq', 'ensembl_{0}.txt'.format(i)).put(Body=b'Test ensembl here!',
                                                                Metadata={'madefor': 'testing', 'quality': quality})

        s3.Object('parseq', 'projects/ensembl.txt').put(Body=b'Test ensembl here!',
                                                    Metadata={'madefor': 'testing', 'quality': 'good'})
        # Setting 'quality' metadata field of files in subdir projections and adding files to mock resource
        for i in range(1,5):
            if i == 1:
                quality = 'good'
            else:
                quality = 'bad'
            s3.Object('parseq', 'projects/ensembl_{0}.txt'.format(i)).put(Body=b'Test ensembl here!',
                                                                Metadata={'madefor': 'testing', 'quality': quality})

        projection_configuration = PrototypeDeserializer('tests/test_s3.yaml')

        s3_driver = S3Driver('test_id', 'test_key', 'us-west-2', projection_configuration.root_projection_uri)
        cls.s3_projector = Projector(s3_driver, projection_configuration.root_projection_uri,
                                     projection_configuration.prototype_tree)

    @classmethod
    def tearDownClass(cls):
        cls.mock.stop()

    def test_create_projections(self):
        """
        Tests if S3 projector creates projections
        """
        created_projections = [n.get_path() for n in self.s3_projector.projection_tree.get_tree_nodes()]

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(6, len(created_projections),
                         msg='Checking if S3 projector created {} projections.'.format(len(created_projections)))

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