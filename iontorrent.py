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

from projections import Projection, ProjectionDriver, ProjectionTree, Projector, ProjectionPrototype, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE
from tests.mock import TorrentSuiteMock

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
        self.authenticate(self.host_url, user, password)

    def __prepare_uri(self, uri):
        """
        Adds appropriate prefix to uri according to uri context
        :param uri: URI string
        :return: URI string
        """
        if re.match('/rundb/api/v1/', uri):
            return uri.replace('/rundb/api/v1/', self.api_url)
        elif re.match('/auth/output/Home/', uri):
            return uri.replace('/auth/output/Home/', self.files_url)
        else:
            return self.api_url + uri

    def authenticate(self, host_url, user, password):
        """
        Creates authorization handler for driver.
        :param host_url: URL of host string
        :param user: user name string
        :param password: password string
        """
        password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        password_manager.add_password(None, host_url, user, password)

        handler = urllib.request.HTTPBasicAuthHandler(password_manager)

        # create "opener" (OpenerDirector instance)
        opener = urllib.request.build_opener(handler)

        # use the opener to fetch a URL
        opener.open(host_url)

        # Install the opener.
        # Now all calls to urllib.request.urlopen use our opener.
        urllib.request.install_opener(opener)

    def get_content(self, uri):
        """
        Opens URI and its contents
        :param uri: URI string
        :return: dict of URI contents
        """

        uri = self.__prepare_uri(uri)

        with urllib.request.urlopen(uri) as f:
            if re.search('\.bam$', uri) or re.search('\.vcf$', uri) or re.search('\.bed$', uri):
                return b''
            else:
                return json.loads(f.readall().decode('utf-8'))

    def load_content(self, uri):
        """
        Load uri contents
        :param uri: URI string
        :return: content bytes
        """
        uri = self.__prepare_uri(uri)
        with urllib.request.urlopen(uri) as f:
            return f.readall()


class TorrentSuiteProjector(Projector):
    def __init__(self, driver, prototype_tree):
        """
        Initializes Torrent Suite Projector with driver, assigns root projection, builds prototype and projection tree.
        :param driver: instance of TorrentSuiteDriver
        :param prototype_tree: tree of ProjectionPrototype objects to build projection upon
        """
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree()
        self.root_projection = Projection('/', 'experiment?status=run&limit=1&order_by=-id')
        self.projection_tree.add_projection(self.root_projection, None)

        self.create_projection_tree({'/': prototype_tree},
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

        content = self.driver.load_content(uri)
        logger.info('Got path content: %s\n', path)

        self.projections[path].size = len(content)

        file_header = 3
        resource_io = io.BytesIO(content)

        return file_header, resource_io


# For smoke testing
def main(mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)
    mock_torrent_suite = TorrentSuiteMock('mockiontorrent.com', 'tests/mock_resource')

    projection_configuration = PrototypeDeserializer('torrent_suite_config.yaml')

    projection_dirver = TorrentSuiteDriver(projection_configuration.resource_uri, 'ionadmin', '0ECu1lW')
    projection_filesystem.projection_manager = TorrentSuiteProjector(projection_dirver,
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



