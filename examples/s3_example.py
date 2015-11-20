import argparse
from aws_s3 import S3Driver
import logging
import logging.config
from fuse import FUSE
from projections import PrototypeDeserializer, Projector
from moto import mock_s3
import boto3
from filesystem import ProjectionFilesystem


def main(cfg_path, mountpoint, data_folder, foreground=True):

    # Setting up mock S3 resource
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

    # Then mock S3 resource is prepared we initialize projection filesystem
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)

    # PrototypeDeserializer class builds prototype tree using which projections will be created,
    # specifies root projection uri and resource uri which will be projected
    projection_configuration = PrototypeDeserializer(cfg_path)

    # This driver is used to connect with projection resource, in case of S3 we specify id and .
    projection_driver = S3Driver('test_id', 'test_key',
                                 'us-west-2', projection_configuration.root_projection_uri)

    projection_filesystem.projection_manager = Projector(projection_driver,
                                                         projection_configuration.root_projection_uri,
                                                         projection_configuration.prototype_tree).projection_tree

    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='S3 projection example.')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory)


