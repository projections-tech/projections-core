import logging
import logging.config
import time
from unittest import TestCase, skip
from projections import PrototypeDeserializer, Projection
import aws_s3

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('s3_test')

class TestS3Driver(TestCase):
    def setUp(self):
        self.driver = aws_s3.S3Driver('parseq')

    def test_bucket_contents(self):
        """
        Tests if driver gets proper bucket on init
        """
        exp_contents = ['projects/', 'projects/ensembl.txt']
        for obj in exp_contents:
            self.assertIn(obj, [o.key for o in self.driver.bucket.objects.all()], msg='Checking objects existence.')

    def test_load_uri_contents_stream(self):
        """
        Tests driver load_uri_contents_stream method
        """
        with open('tests/test_ensembl.txt', 'rb') as t_f:
            test_ensembl = t_f.read()
        for key, exp_content in {'projects/': b'', 'projects/ensembl.txt': test_ensembl}.items():
            self.assertEqual(exp_content, self.driver.get_uri_contents_stream(key),
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
