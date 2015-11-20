#!/usr/bin/env python3

import logging
import logging.config
import os
import boto3
import argparse

from projections import ProjectionDriver, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE
from moto import mock_s3

logger = logging.getLogger('s3_projection')


class S3Driver(ProjectionDriver):

    def __init__(self, aws_access_key_id, aws_secret_access_key, region_name, bucket_name):
        """
        Initialize driver which will be used to interact with host.
        :param bucket_name: name of projected S3 bucket
        """

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

# For smoke testing
def main(mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    mock = mock_s3()
    mock.start()
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

    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)

    projection_configuration = PrototypeDeserializer('aws_s3.yaml')

    sra_driver = S3Driver('test_id', 'test_key',
                          'us-west-2', projection_configuration.root_projection_uri)
    projection_filesystem.projection_manager = Projector(sra_driver, projection_configuration.root_projection_uri,
                                                         projection_configuration.prototype_tree).projection_tree
    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    parser = argparse.ArgumentParser(description='AWS S3 projection.')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    args = parser.parse_args()

    main(args.mount_point, args.data_directory)
