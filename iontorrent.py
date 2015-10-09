#!/usr/bin/env python3
__author__ = 'abragin'

import logging
import logging.config
import io
import re
import json
import os
import stat
import sys
import time
import urllib.request
from urllib.parse import urljoin

from projections import Projection,  ProjectionDriver, ProjectionTree, Projector, ProjectionPrototype
from filesystem import ProjectionFilesystem
from fuse import FUSE
from tests.mock import TorrentSuiteMock

logger = logging.getLogger('iontorrent_projection')


class TorrentSuiteDriver(ProjectionDriver):
    def __init__(self, host_url, user, password):
        self.host_url = 'http://{}'.format(host_url)
        self.api_url = 'http://{}/rundb/api/v1/'.format(host_url)
        self.files_url = urljoin(self.host_url, '/auth/output/Home/')
        self.authenticate(self.host_url, user, password)

    def authenticate(self, host_url, user, password):
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
        with urllib.request.urlopen(uri) as f:
            if re.search('\.bam$', uri) or re.search('\.vcf$', uri) or re.search('\.bed$', uri):
                return b''
            else:
                return json.loads(f.readall().decode('utf-8'))
        # TODO write this function in a way which automatically assigns right host address according to request

    def load_content(self, uri):
        with urllib.request.urlopen(uri) as f:
            return f.readall()

class TorrentSuiteProjector(Projector):
    def __init__(self, driver):
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree()
        self.root_projection = Projection('/', self.driver.api_url + 'experiment?status=run&limit=1&order_by=-id')
        self.projection_tree.add_projection(self.root_projection, None)

        prototypes = self.prepare_prototypes()
        self.create_projection_tree(prototypes, projection_tree=self.projection_tree, parent_projection=self.root_projection)
        self.projections = self.projection_tree.projections

    def prepare_prototypes(self):
        """
        This method is STUB, in future prototype deserialization will be implemented
        :return: root prototype object
        """
        # TODO implement resource uri resolution in driver
        experiment_prototype = ProjectionPrototype('directory')
        experiment_prototype.name = "content['displayName'].replace(' ', '_')"
        experiment_prototype.uri = '["{0}"+object["resource_uri"] for object in environment["objects"]]'.format(self.driver.host_url)

        exp_metadata_prototype = ProjectionPrototype('file', experiment_prototype)
        exp_metadata_prototype.name = "'metadata.json'"
        exp_metadata_prototype.uri = '["{0}" + environment["resource_uri"]]'.format(self.driver.host_url)

        exp_plan_metadata_prototype = ProjectionPrototype('file', experiment_prototype)
        exp_plan_metadata_prototype.name = "'plannedexperiment.json'"
        exp_plan_metadata_prototype.uri = '["{0}" + environment["plan"]]'.format(self.driver.host_url)

        result_prototype = ProjectionPrototype('directory', experiment_prototype)
        result_prototype.name = "path.split(content['filesystempath'])[1]"
        result_prototype.uri = "['{0}' + res for res in environment['results']]".format(self.driver.host_url)

        sample_prototype = ProjectionPrototype('directory', result_prototype)
        sample_prototype.name = "content['name']"
        sample_prototype.uri = "['{0}' + sample['resource_uri'] for sample in context[1]['samples']]".format(self.driver.host_url)

        sample_metadata_prototype = ProjectionPrototype('file', experiment_prototype)
        sample_metadata_prototype.name = "'metadata.json'"
        sample_metadata_prototype.uri = '["{0}" + environment["resource_uri"]]'.format(self.driver.host_url)

        bam_prototype = ProjectionPrototype('file', sample_prototype)
        bam_prototype.name = "environment['name']+'.bam'"
        bam_prototype.uri = "['{0}'+ path.split(context[2]['filesystempath'])[1] + '/'" \
                            "+fetch_context('{1}'+context[1]['plan'])['barcodedSamples'][environment['name']]['barcodes'][0]" \
                            " + '_rawlib.bam']".format(self.driver.files_url, self.driver.host_url)

        plugin_result_prototype = ProjectionPrototype('directory', sample_prototype)
        plugin_result_prototype.name = "path.basename(content['path'])"
        plugin_result_prototype.uri = "['{0}' + p_res for p_res in context[2]['pluginresults']" \
                                      " if 'variant' in fetch_context('{0}' + p_res)['pluginName'] and 'VFNA' not in fetch_context('{0}' + p_res)['pluginName']]".format(self.driver.host_url)

        bed_prototype = ProjectionPrototype('file', plugin_result_prototype)
        bed_prototype.name = "path.basename(environment['store']['targets_bed'])"
        bed_prototype.uri = "[ '{0}' + path.basename(context[2]['filesystempath']) + '/plugin_out/'" \
                                " + path.basename(environment['path']) " \
                                " + '/' + path.basename(environment['store']['targets_bed'])]".format(self.driver.files_url,
                                                                                                self.driver.host_url)

        local_settings_prototype = ProjectionPrototype('file', plugin_result_prototype)
        local_settings_prototype.name = " 'variant_caller_settings.json' "
        local_settings_prototype.uri = "[ '{0}' + path.basename(context[2]['filesystempath']) + '/plugin_out/'" \
                                " + path.basename(environment['path']) " \
                                " + '/' + 'local_parameters.json' ]".format(self.driver.files_url, self.driver.host_url)

        for variant_file_name in ['TSVC_variants.vcf', 'all.merged.vcf', 'indel_assembly.vcf',
                                                      'indel_variants.vcf', 'small_variants.left.vcf',
                                                      'small_variants.vcf', 'small_variants_filtered.vcf',
                                                      'small_variants.sorted.vcf', 'SNP_variants.vcf']:
            vcf_prototype = ProjectionPrototype('file', plugin_result_prototype)
            vcf_prototype.name = "'{0}'".format(variant_file_name)
            vcf_prototype.uri = "['{0}' + path.basename(context[2]['filesystempath']) + '/plugin_out/'" \
                                " + path.basename(environment['path']) " \
                                " + '/' + fetch_context('{1}'" \
                                " + context[1]['plan'])['barcodedSamples'][context[3]['name']]['barcodes'][0] " \
                                " + '/{2}']".format(self.driver.files_url,
                                                    self.driver.host_url,
                                                    variant_file_name)
            plugin_result_prototype.children[vcf_prototype.name] = vcf_prototype

        experiment_prototype.children[result_prototype.name] = result_prototype
        experiment_prototype.children[exp_metadata_prototype.name] = exp_metadata_prototype
        experiment_prototype.children[exp_plan_metadata_prototype] = exp_plan_metadata_prototype

        result_prototype.children[sample_prototype.name] = sample_prototype

        sample_prototype.children[sample_metadata_prototype.name] = sample_metadata_prototype
        sample_prototype.children[bam_prototype.name] = bam_prototype
        sample_prototype.children[plugin_result_prototype.name] = plugin_result_prototype

        plugin_result_prototype.children[local_settings_prototype.name] = local_settings_prototype
        plugin_result_prototype.children[bed_prototype.name] = bed_prototype

        return {'/': experiment_prototype}

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
    mock_url = mock_torrent_suite.mock_url

    projection_dirver = TorrentSuiteDriver(mock_url, 'ionadmin', '0ECu1lW')
    projection_filesystem.projection_manager = TorrentSuiteProjector(projection_dirver)
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




