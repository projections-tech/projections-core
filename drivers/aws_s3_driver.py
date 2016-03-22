#!/usr/bin/env python3

import logging
import logging.config

import boto3

from projections import ProjectionDriver

logger = logging.getLogger('s3_driver')


class S3Driver(ProjectionDriver):
    def __init__(self, bucket_name, driver_config_path, script_dir):
        """
        Initialize driver which will be used to interact with host.
        :param bucket_name: name of projected S3 bucket
        """
        self.driver_configuration = self.read_config(script_dir, driver_config_path)

        aws_access_key_id, aws_secret_access_key, region_name = self.driver_configuration['aws_access_key_id'], \
                                                                self.driver_configuration['aws_secret_access_key'], \
                                                                self.driver_configuration['region_name']

        session = boto3.session.Session(aws_access_key_id=aws_access_key_id,
                                        aws_secret_access_key=aws_secret_access_key,
                                        region_name=region_name)
        self.s3_resource = session.resource('s3')
        self.bucket_name = bucket_name
        self.bucket = self.s3_resource.Bucket(bucket_name)
        logger.debug('Current bucket: %s', self.bucket)

    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """
        metadata = dict()

        if uri != self.bucket_name:
            current_object = self.bucket.Object(key=uri)
            metadata['name'] = current_object.key
            metadata['metadata'] = current_object.metadata
            metadata['size'] = current_object.content_length
            metadata['content_encoding'] = current_object.content_encoding
            metadata['content_type'] = current_object.content_type
            metadata['resource_uri'] = uri
            return metadata
        else:
            return {'bucket_contents': [self.get_uri_contents_as_dict(o.key) for o in self.bucket.objects.all()]}

    def get_uri_contents_as_bytes(self, uri):
        """
        Load uri contents as bytes
        :param uri: URI string
        :return: content bytes
        """
        return self.bucket.Object(key=uri).get()['Body'].read()