#!/usr/bin/env python3
__author__ = 'abragin'

import argparse
import logging
import logging.config
import re
import json
import os
import urllib.request
from urllib.parse import urljoin

from projections import ProjectionDriver, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE
from tests.mock import MockResource

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
                return {}
            else:
                meta = json.loads(f.readall().decode('utf-8'))

                # Here we add required for projection fields to metadata
                if 'resource_uri' in meta:
                    if '/experiment/' in meta['resource_uri']:
                        # Experiment is hierarchically upper element, it`s fields 'plan' and 'results' used by lower
                        # level projections.
                        # Firstly we add 'plan' field to metadata
                        plan_uri = self.__prepare_uri(meta['plan'])
                        with urllib.request.urlopen(plan_uri) as p_f:
                            meta['plan'] = json.loads(p_f.readall().decode('utf-8'))
                        # Secondly 'results' field
                        temp = []
                        for res in meta['results']:
                            res_uri = self.__prepare_uri(res)
                            with urllib.request.urlopen(res_uri) as r:
                                temp.append(json.loads(r.readall().decode('utf-8')))
                        meta['results'] = temp

                    elif '/results/' in meta['resource_uri']:
                        # Here we resolve 'experiment' of upper level resource using recursive call
                        exp_uri = self.__prepare_uri(meta['experiment'])
                        meta['experiment'] = self.get_uri_contents_as_dict(exp_uri)
                        # Adding 'pluginresults' to result meta
                        temp = []
                        for p_res in meta['pluginresults']:
                            p_res_uri = self.__prepare_uri(p_res)
                            with urllib.request.urlopen(p_res_uri) as p_r:
                                temp.append(json.loads(p_r.readall().decode('utf-8')))
                        meta['pluginresults'] = temp
                return meta

    def get_uri_contents_as_bytes(self, uri):
        """
        Get uri contents as bytes
        :param uri: URI string
        :return: content bytes
        """
        uri = self.__prepare_uri(uri)
        with urllib.request.urlopen(uri) as f:
            return f.readall()


# For smoke testing
def main(cfg_path, mountpoint, data_folder, projection_name, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)
    mock_torrent_suite = MockResource('tests/torrent_suite_mock.json')
    projection_configuration = PrototypeDeserializer(cfg_path)
    projection_driver = TorrentSuiteDriver(projection_configuration.resource_uri, 'ionadmin', '0ECu1lW')

    from db_projector import DBProjector
    ion_torrent_projection_tree = DBProjector(projection_name, projection_driver, 'viktor', 'test',
                                              projection_configuration.prototype_tree,
                                              projection_configuration.root_projection_uri)
    projection_filesystem.projection_manager = ion_torrent_projection_tree
    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    parser = argparse.ArgumentParser(description='Torrent Suite projection.')
    parser.add_argument('-p', '--projection-name', required=True, help='name of current projection')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory, args.projection_name)




