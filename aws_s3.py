#!/usr/bin/env python3

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

from projections import Projection, ProjectionDriver, ProjectionTree, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE

logger = logging.getLogger('s3_projection')

class S3Driver(ProjectionDriver):
    def __init__(self, bucket_name):
        """
        Initialize driver which will be used to interact with host.
        :param bucket_name: name of projected S3 bucket
        """
        session = boto3.session.Session(aws_access_key_id='AKIAIONUXTO6TR3UU3TQ',
                                        aws_secret_access_key='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',
                                        region_name='us-east-1')
        self.s3_resource = session.resource('s3')
        self.bucket = self.s3_resource.Bucket(bucket_name)
        logger.debug('Bucket initialized: %s', self.bucket)

    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """
        current_object = self.bucket.Object(key=uri)
        metadata = dict()

        metadata['metadata'] = current_object.metadata
        metadata['size'] = current_object.content_length
        metadata['content_encoding'] = current_object.content_encoding
        metadata['content_type'] = current_object.content_type
        return metadata

    def get_uri_contents_stream(self, uri):
        """
        Load uri contents
        :param uri: URI string
        :return: content bytes
        """
        return self.bucket.Object(key=uri).get()['Body'].read()



