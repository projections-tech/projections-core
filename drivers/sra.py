#!/usr/bin/env python3

import argparse
import json
import logging
import logging.config
import os
import re
import subprocess

import xmltodict
from Bio import Entrez

from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import ProjectionDriver, Projector, PrototypeDeserializer
from tests.mock import MockResource

logger = logging.getLogger('sra_projection')


class SRADriver(ProjectionDriver):
    def __init__(self, uri, driver_config, script_dir):
        """
        Initialize driver telling NCBI using user email and name of program
        :param driver_config: driver configuration dictionary
        :return:
        """
        Entrez.email = driver_config['email']
        Entrez.tool = 'sra_projection_manager'
        self.driver_cache = {}

    def get_uri_contents_as_dict(self, uri):
        """
        Loads content from SRA using Biopython in driver cache
        :param uri: str containing uri
        :return: dict of query contents
        """
        # Info about Biopython`s eutils: http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter:entrez
        # Query looks as: 'query:Test_species'
        logger.debug('Current query: %s', uri)

        if uri not in self.driver_cache:
            uri_parts = uri.split(':')
            # Returns esearch response dict for SRA database.
            esearch_handle = Entrez.esearch(db='sra', term=uri_parts[1], retmax=uri_parts[2])
            search_result = Entrez.read(esearch_handle)
            # Id`s of experiment resources
            for sra_id in search_result['IdList']:
                # Driver accesses experiments using id`s like sra_id:1214564
                res_uri = 'sra_id:{0}'.format(sra_id)
                # Handler to fetch SRA database by experiment id
                fetch_handler = Entrez.efetch(db='sra', id=sra_id)

                # Converts nested ordered dicts to default dicts required by ObjectPath
                experiment = json.loads(json.dumps(xmltodict.parse(fetch_handler.read())))
                # Setting resource uri for experiment on driver
                experiment['resource_uri'] = res_uri

                # Experiment resource contains data about runs, which contain id`s of SAM files in database
                # Run set is most times dict, but sometimes list, treating dict as list to resolve inconsistency
                run_set = experiment['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN']
                if not isinstance(run_set, list):
                    run_set = [run_set]
                # Setting corrected run set field for experiment
                experiment['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN'] = run_set
                # Adding experiment resource in driver cache by uri
                self.driver_cache[res_uri] = experiment

                # Adding runs to driver cache by their accession
                for run in run_set:
                    self.driver_cache[run['@accession']] = run

            # Adding sra_id prefix by which experiment resources will be accessed using driver cache
            search_result['IdList'] = [''.join(['sra_id:', id]) for id in search_result['IdList']]
            # Adding resource uri to search results meta
            search_result['resource_uri'] = uri
            # Converting search result into proper dict, not biopython subclass, because ObjectPath can handle only
            # standard dicts
            search_result = dict(search_result)
            # Adding results of search to driver cache
            self.driver_cache[uri] = search_result
            return search_result
        else:
            return self.driver_cache[uri]

    def load_uri_contents_as_bytes(self, uri):
        """
        Returns stream of uri contents
        :param uri: uri of resource
        :return: stream of resource contents
        """
        logger.debug('Loading query: %s', uri)
        # Regex matches id`s in SRA database
        sam_uri_regex = '(SRR|SRX|ERX|DRX|DRR|ERR)\d+'
        if re.match(sam_uri_regex, uri):
            # Using subprocess.check to run sam-dump, which returns stream after loading of sam file,
            # this approach may be slow, other approaches lock script execution, need to reconsider
            return subprocess.check_output(['./sratoolkit.2.5.4-1-ubuntu64/bin/sam-dump', uri])
        else:
            return json.dumps(self.driver_cache[uri]).encode()


# For smoke testing
def main(cfg_path, mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)
    mock_resource = MockResource('tests/sra_mock.json')

    projection_configuration = PrototypeDeserializer(cfg_path)

    sra_driver = SRADriver('vsvekolkin@parseq.pro')

    sra_projection_tree = Projector(sra_driver, projection_configuration.root_projection_uri,
                                    projection_configuration.prototype_tree).projection_tree

    projection_filesystem.projection_manager = sra_projection_tree
    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    parser = argparse.ArgumentParser(description='SRA projection.')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory)
