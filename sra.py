#!/usr/bin/env python3

import logging
import logging.config
import json
import os
import xmltodict
import subprocess
import argparse
from Bio import Entrez
from tests.sra_mock import SRAMock


from projections import ProjectionDriver, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE


logger = logging.getLogger('sra_projection')


class SRADriver(ProjectionDriver):
    def __init__(self, email):
        Entrez.email = email
        Entrez.tool = 'sra_projection_manager'
        self.query_cache = {}

    def get_uri_contents_as_dict(self, query):
        """
        Loads content from SRA using Biopython in driver cache, with queries in format: "query_type:query"
        :param query: str containing query to SRA
        :return: dict of query contents
        """
        query = query.split(':')
        logger.debug('Current query: %s', query)
        query_type = query[0]
        # Info about Biopython`s eutils: http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter:entrez
        # Query looks as: 'query:Test_species'
        if query_type == 'query':
            # Returns esearch response dict.
            esearch_handle = Entrez.esearch(db='sra', term=query[1], retmax=query[2])
            return Entrez.read(esearch_handle)
        # Query looks as 'search_id:102354'
        elif query_type == 'search_id':
            if query[1] not in self.query_cache:
                fetch_handler = Entrez.efetch(db='sra', id=query[1])
                sample_dict = xmltodict.parse(fetch_handler.read())
                search_query_contents = sample_dict['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
                # Adding resources to cache dict with their IDs as keys
                self.query_cache[query[1]] = search_query_contents
                self.query_cache[search_query_contents['EXPERIMENT']['@accession']] = search_query_contents
                # Run set is most times dict, but sometimes list, treating dict as list to resolve inconsistency
                if isinstance(search_query_contents['RUN_SET']['RUN'], list):
                    for run in search_query_contents['RUN_SET']['RUN']:
                        self.query_cache[run['@accession']] = run['@accession']
                else:
                    self.query_cache[search_query_contents['RUN_SET']['RUN']['@accession']] = search_query_contents['RUN_SET']['RUN']['@accession']
                return search_query_contents
            else:
                return self.query_cache[query[1]]
        elif query_type == 'get_experiment_runs':
            experiment_data = self.query_cache[query[1]]
            # Run set is most times dict, but sometimes list, treating dict as list to resolve inconsistency
            if isinstance(experiment_data['RUN_SET']['RUN'], list):
                return [run['@accession'] for run in experiment_data['RUN_SET']['RUN']]
            else:
                return [experiment_data['RUN_SET']['RUN']['@accession']]
        else:
            return self.query_cache[query[1]]

    def get_uri_contents_as_bytes(self, query):
        """
        Open URI and return bytes massive
        :param uri: URI string
        :return: bytes massive
        """
        logger.debug('Loading query: %s', query)
        query = query.split(':')
        query_type = query[0]
        if query_type == 'load_run':
            return subprocess.check_output(['./sratoolkit.2.5.4-1-ubuntu64/bin/sam-dump', query[1]])
        else:
            return json.dumps(self.query_cache[query[1]]).encode()


# For smoke testing
def main(cfg_path, mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)
    mock_resource = SRAMock('http://eutils.ncbi.nlm.nih.gov', 'tests/mock_resource')

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
