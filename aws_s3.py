#!/usr/bin/env python3

import logging
import logging.config
import io
import os
import time
import boto3
import argparse

from projections import Projection, ProjectionDriver, ProjectionTree, Projector, PrototypeDeserializer
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
            return metadata
        else:
            return [o.key for o in self.bucket.objects.all()]

    def get_uri_contents_as_stream(self, uri):
        """
        Load uri contents
        :param uri: URI string
        :return: content bytes
        """
        return self.bucket.Object(key=uri).get()['Body'].read()


class S3Projector(Projector):
    def __init__(self, driver, root_projection, root_prototype):
        """
        Initializes S3 Projector with driver, assigns root projection, builds prototype and projection tree.
        :param driver: instance of S3Driver
        :param prototype_tree: tree of ProjectionPrototype objects to build projection upon
        """
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree()
        self.root_projection = root_projection
        self.projection_tree.add_projection(self.root_projection, None)

        self.create_projection_tree({'/': root_prototype},
                                    projection_tree=self.projection_tree,
                                    parent_projection=self.root_projection)
        self.projections = self.projection_tree.projections

    def is_managing_path(self, path):
        return path in self.projections

    def get_projections(self, path):
        logger.info('Requesting projections for path: %s', path)
        projections = []
        for p in self.projections:
            if p.startswith(path):
                # limit projections to one level only
                logging.debug('Analyzing projection candidate: %s', p)
                suffix = p[len(path):]
                if suffix and suffix[0] == '/':
                    suffix = suffix[1:]
                logger.debug('Path suffix: %s', suffix)
                if suffix and '/' not in suffix:
                    projections.append('/' + suffix)

        logger.debug('Returning projections: %s', projections)
        return projections

    def get_attributes(self, path):
        assert path in self.projections

        projection = self.projections[path]

        now = time.time()
        attributes = dict()

        # Set projection attributes

        # This is implementation specific and should be binded to projector data
        attributes['st_atime'] = now
        # This may be implemented as last projection cashing time is casing is enabled
        attributes['st_mtime'] = now
        # On Unix this is time for metedata modification we can use the same conception
        attributes['st_ctime'] = now
        # If this is projection the size is zero
        attributes['st_size'] = projection.size
        # Set type to link anf grant full access to everyone
        attributes['st_mode'] = (projection.type | 0o0777)
        # Set number of hard links to 0
        attributes['st_nlink'] = 0
        # Set id as inode number.
        attributes['st_ino'] = 1

        return attributes

    def open_resource(self, path):
        uri = self.projections[path].uri

        content = self.driver.get_uri_contents_as_stream(uri)
        logger.info('Got path content: %s\n', path)

        self.projections[path].size = len(content)

        file_header = 3
        resource_io = io.BytesIO(content)

        return file_header, resource_io


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

    root_projection = Projection('/', projection_configuration.root_projection_uri)

    sra_driver = S3Driver('test_id', 'test_key',
                          'us-west-2', projection_configuration.root_projection_uri)
    projection_filesystem.projection_manager = S3Projector(sra_driver, root_projection,
                                                           projection_configuration.prototype_tree)
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
