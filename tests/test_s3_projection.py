import logging
import logging.config
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
        exp_contents = ['projects/', 'projects/ensembl.txt']
        for obj in exp_contents:
            self.assertIn(obj, self.driver.bucket_contents.keys())

    def test_load_uri_contents_stream(self):
        with open('tests/test_ensembl.txt', 'rb') as t_f:
            test_ensembl = t_f.read()
        for key, exp_content in {'projects/':b'', 'projects/ensembl.txt':test_ensembl}.items():
            self.assertEqual(exp_content, self.driver.get_uri_contents_stream(key))

    def test_get_uri_contents_as_dict(self):
        exp_contents = ['projects/', 'projects/ensembl.txt']
        for key in exp_contents:
            logger.debug(self.driver.get_uri_contents_as_dict(key))
