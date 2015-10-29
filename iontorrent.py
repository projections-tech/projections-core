#!/usr/bin/env python3
__author__ = 'abragin'

import logging
import logging.config
import re
import json
import os
import sys
import urllib.request
from urllib.parse import urljoin

from projections import ProjectionDriver, Projector, PrototypeDeserializer
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
        if re.match('/rundb/api/v1/', uri):
            return uri.replace('/rundb/api/v1/', self.api_url)
        elif re.match('/auth/output/Home/', uri):
            return uri.replace('/auth/output/Home/', self.files_url)
        else:
            return self.api_url + uri

    def authenticate(self, user, password):
        """
        Creates authorization handler for driver.
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
        Open URI and return dict of its contents
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

    def get_uri_contents_as_stream(self, uri):
        """
        Get uri contents as stream
        :param uri: URI string
        :return: content bytes
        """
        uri = self.__prepare_uri(uri)
        with urllib.request.urlopen(uri) as f:
            return f.readall()


# For smoke testing
def main(mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)
    mock_torrent_suite = TorrentSuiteMock('mockiontorrent.com', 'tests/mock_resource/torrent_suite_mock_data')

    projection_configuration = PrototypeDeserializer('torrent_suite_config.yaml')
    projection_driver = TorrentSuiteDriver(projection_configuration.resource_uri, 'ionadmin', '0ECu1lW')
    projection_filesystem.projection_manager = Projector(projection_driver,
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



