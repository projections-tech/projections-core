import logging
import logging.config
import time
from unittest import TestCase, skip
from projections import PrototypeDeserializer, Projection
from aws_s3 import S3Driver, S3Projector

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('s3_test')

class TestS3Driver(TestCase):
    def setUp(self):
        aws_access_key_id='AKIAIONUXTO6TR3UU3TQ'
        aws_secret_access_key='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'
        region_name='us-east-1'
        self.driver = S3Driver(aws_access_key_id, aws_secret_access_key, region_name, 'parseq')

    def test_bucket_contents(self):
        """
        Tests if driver gets proper bucket on init
        """
        exp_contents = ['projects/', 'projects/ensembl.txt']
        for obj in exp_contents:
            self.assertIn(obj, [o.key for o in self.driver.bucket.objects.all()], msg='Checking objects existence.')

    def test_load_uri_contents_as_stream(self):
        """
        Tests driver load_uri_contents_stream method
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
        exp_contents = ['projects/', 'projects/ensembl.txt']
        contents_meta = [{'content_encoding': None, 'size': 0,
                          'content_type': 'binary/octet-stream', 'metadata': {}},
                         {'content_encoding': None, 'size': 1387,
                          'content_type': 'text/plain', 'metadata': {'madefor': 'testing'}}]
        for key, meta in zip(exp_contents, contents_meta):
            self.assertDictEqual(meta, self.driver.get_uri_contents_as_dict(key),
                                 msg='Checking meta of object with URI: {0}'.format(key))
