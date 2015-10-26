import logging
import logging.config
import io
import re
import json
import os
import sys
import time
import boto3
import urllib.request
from urllib.parse import urljoin

from projections import Projection, ProjectionDriver, ProjectionTree, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE

logger = logging.getLogger('s3_projection')

class S3Driver(ProjectionDriver):
    def __init__(self, bucket_names):
        """
        Initialize driver which will be used to interact with host.
        :param host_url: URL of host string
        """
        session = boto3.session.Session(aws_access_key_id='AKIAIONUXTO6TR3UU3TQ',
                                        aws_secret_access_key='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',
                                        region_name='us-east-1')
        self.s3_resource = session.resource('s3')


    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """

        with urllib.request.urlopen(uri) as f:
            if re.search('\.bam$', uri) or re.search('\.vcf$', uri) or re.search('\.bed$', uri):
                return b''
            else:
                return json.loads(f.readall().decode('utf-8'))

    def load_uri_contents_stream(self, uri):
        """
        Load uri contents
        :param uri: URI string
        :return: content bytes
        """
        uri = self.__prepare_uri(uri)
        with urllib.request.urlopen(uri) as f:
            return f.readall()

if __name__ == '__main__':
    test_driver = S3Driver(['parseq'])
