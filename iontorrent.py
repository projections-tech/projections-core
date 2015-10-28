#!/usr/bin/env python3
__author__ = 'abragin'

import logging
import logging.config
import io
import re
import json
import os
import sys
import time
import urllib.request
from urllib.parse import urljoin

from projections import Projection, ProjectionDriver, ProjectionTree, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE
from tests.torrent_suite_mock import TorrentSuiteMock

logger = logging.getLogger('iontorrent_projection')


class TorrentSuiteDriver(ProjectionDriver):
    def __init__(self, host_url, user, password):
        """
        Initialize driver which will be used to interact with host.
        :param host_url: URL of host string
        :param user: user name string
        :param password: password string
        """
        self.host_url = 'http://{}'.format(host_url)
        self.api_url = 'http://{}/rundb/api/v1/'.format(host_url)
        self.files_url = urljoin(self.host_url, '/auth/output/Home/')
        self.authenticate(user, password)

    def __prepare_uri(self, uri):
        """
        Adds appropriate prefix to uri according to uri context
        :param uri: URI string
        :return: URI string
        """
        logger.debug('Driver URI: %s', type(uri))
        if re.match('/rundb/api/v1/', uri):
            return uri.replace('/rundb/api/v1/', self.api_url)
        elif re.match('/auth/output/Home/', uri):
            return uri.replace('/auth/output/Home/', self.files_url)
        else:
            return self.api_url + uri

    def authenticate(self, user, password):
        """
        Creates authorization handler for driver.
        :param host_url: URL of host string
        :param user: user name string
        :param password: password string
        """
        password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        password_manager.add_password(None, self.host_url, user, password)

        handler = urllib.request.HTTPBasicAuthHandler(password_manager)

        # create "opener" (OpenerDirector instance)
        opener = urllib.request.build_opener(handler)

        # use the opener to fetch a URL
        opener.open(self.host_url)

        # Install the opener.
        # Now all calls to urllib.request.urlopen use our opener.
        urllib.request.install_opener(opener)

    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """

        uri = self.__prepare_uri(uri)
        # TODO move this functionality to caller function
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


class TorrentSuiteProjector(Projector):
    def __init__(self, driver, root_projection_uri, prototype_tree):
        """
        Initializes Torrent Suite Projector with driver, assigns root projection, builds prototype and projection tree.
        :param driver: instance of TorrentSuiteDriver
        :param prototype_tree: tree of ProjectionPrototype objects to build projection upon
        """
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree(p_name='/', p_uri=root_projection_uri)

        self.create_projection_tree({'/': prototype_tree},
                                    projection_tree=self.projection_tree)

    def is_managing_path(self, path):
        if self.projection_tree.get_projection(path):
            return True
        else:
            return False

    def get_projections(self, path):
        logger.info('Requesting projections for path: %s', path)
        projections = [c.projection.path for c in self.projection_tree.get_children(path)]
        logger.info('Returning projections: %s', projections)
        return projections

    def get_attributes(self, path):
        assert self.projection_tree.get_projection(path) is not None

        projection = self.projection_tree.get_projection(path).projection

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
        projection_on_path = self.projection_tree.get_projection(path).projection
        uri = projection_on_path.uri

        content = self.driver.load_uri_contents_stream(uri)
        logger.info('Got path content: %s\n', path)

        projection_on_path.size = len(content)

        file_header = 3
        resource_io = io.BytesIO(content)

        return file_header, resource_io


# For smoke testing
def main(mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)
    mock_torrent_suite = TorrentSuiteMock('mockiontorrent.com', 'tests/mock_resource/torrent_suite_mock_data')

    projection_configuration = PrototypeDeserializer('torrent_suite_config.yaml')
    projection_dirver = TorrentSuiteDriver(projection_configuration.resource_uri, 'ionadmin', '0ECu1lW')
    projection_filesystem.projection_manager = TorrentSuiteProjector(projection_dirver,
                                                                     projection_configuration.root_projection_uri,
                                                                     projection_configuration.prototype_tree)
    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    # TODO: replace with argparse
    # Get mount point from args
    if len(sys.argv) != 3:
        print('usage: %s <mountpoint> <data folder>' % sys.argv[0])
        exit(1)

    main(sys.argv[1], sys.argv[2])



