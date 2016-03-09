#!/usr/bin/env python3

import argparse
import logging
import logging.config
import os

from Bio import Entrez

from filesystem import ProjectionFilesystem
from fuse import FUSE
from projections import ProjectionDriver, PrototypeDeserializer
from tests.mock import MockResource

logger = logging.getLogger('genbank_projection')


class GenbankDriver(ProjectionDriver):
    def __init__(self, email):
        Entrez.email = email
        Entrez.tool = 'genbank_projection_manager'
        self.driver_cache = {}

    def get_uri_contents_as_dict(self, uri):
        """
        Loads content from SRA using Biopython, with queries in format: "query_type:query"
        :param query: str containing query to SRA
        :return: dict of query contents
        """
        if uri not in self.driver_cache:
            query = uri.split(':')
            logger.debug('Current query: %s', query)
            # Info about Biopython`s eutils: http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter:entrez
            # Query looks as: 'query:Test_species'

            # Returns esearch response dict.
            esearch_handle = Entrez.esearch(db='nucleotide', term=query[1], retmax=query[2])
            search_result = dict(Entrez.read(esearch_handle))
            search_result['resource_uri'] = uri
            search_result['fasta_files'] = []
            search_result['gb_files'] = []

            for nuc_id in search_result['IdList']:
                fasta_id = 'fasta:{}'.format(nuc_id)
                gb_id = 'gb:{}'.format(nuc_id)

                search_result['fasta_files'].append(fasta_id)
                search_result['gb_files'].append(gb_id)

                self.driver_cache[fasta_id] = {'resource_uri': fasta_id}
                self.driver_cache[gb_id] = {'resource_uri': gb_id}
            self.driver_cache[uri] = search_result
        return self.driver_cache[uri]

    def get_uri_contents_as_bytes(self, query):
        """
        Open URI and return it`s content as bytes
        :param query: query to driver
        :return: bytes massive
        """
        logger.debug('Loading query: %s', query)
        query = query.split(':')
        query_type = query[0]
        if query_type == 'search_query':
            return Entrez.esearch(db='nuccore', term=query[1], retmax=query[2])
        elif query_type == 'gb':
            # Returns query gb file bytes
            gb = Entrez.efetch(db='nuccore', id=query[1], rettype='gb', retmode='text')
            return gb.read().encode()
        elif query_type == 'fasta':
            # Returns query fasta file bytes
            fasta = Entrez.efetch(db='nuccore', id=query[1], rettype='fasta', retmode='text')
            return fasta.read().encode()


# For smoke testing
def main(config_path, mountpoint, data_folder, projection_name, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)
    projection_configuration = PrototypeDeserializer(config_path)

    genbank_driver = GenbankDriver('vsvekolkin@parseq.pro')

    from db_projector import DBProjector
    projection_filesystem.projection_manager = DBProjector(projection_name, genbank_driver,
                                                           projection_configuration.prototype_tree,
                                                           projection_configuration.root_projection_uri)

    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    mock_resource = MockResource('tests/genbank_mock.json')

    parser = argparse.ArgumentParser(description='Genbank projection.')
    parser.add_argument('-p', '--projection-name', required=True, help='name of current projection')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory, args.projection_name)
    mock_resource.deactivate()
