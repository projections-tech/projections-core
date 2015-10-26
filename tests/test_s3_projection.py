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
