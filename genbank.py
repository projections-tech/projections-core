#!/usr/bin/env python3

import logging
import logging.config
import os
import argparse
from Bio import Entrez
from tests.genbank_mock import GenbankMock

from projections import ProjectionDriver, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE

logger = logging.getLogger('genbank_projection')


class GenbankDriver(ProjectionDriver):
    def __init__(self, email):
        Entrez.email = email
        Entrez.tool = 'genbank_projection_manager'
        self.query_cache = {}

    def get_uri_contents_as_dict(self, query):
        """
        Loads content from SRA using Biopython, with queries in format: "query_type:query"
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
            esearch_handle = Entrez.esearch(db='nucleotide', term=query[1], retmax=query[2])
            return Entrez.read(esearch_handle)
        else:
            return {'id': query[0]}

    def get_uri_contents_as_stream(self, query):
        logger.debug('Loading query: %s', query)
        query = query.split(':')
        query_type = query[0]
        if query_type == 'query':
            return Entrez.esearch(db='nuccore', term=query[1], retmax=query[2])
        elif query_type == 'get_gb':
            # Returns query gb file bytes
            gb = Entrez.efetch(db='nuccore', id=query[1], rettype='gb', retmode='text')
            return gb.read().encode()
        elif query_type == 'get_fasta':
            # Returns query fasta file bytes
            fasta = Entrez.efetch(db='nuccore', id=query[1], rettype='fasta', retmode='text')
            return fasta.read().encode()


# For smoke testing
def main(config_path, mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)
    projection_configuration = PrototypeDeserializer(config_path)

    genbank_driver = GenbankDriver('vsvekolkin@parseq.pro')

    projection_filesystem.projection_manager = Projector(genbank_driver, projection_configuration.root_projection_uri,
                                                         projection_configuration.prototype_tree).projection_tree

    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))
    mock_resource = GenbankMock('http://eutils.ncbi.nlm.nih.gov', 'tests/mock_resource')
    parser = argparse.ArgumentParser(description='Genbank projection.')
    parser.add_argument('-m', '--mount-point', required=True, help='specifies mount point path on host')
    parser.add_argument('-d', '--data-directory', required=True, help='specifies data directory path on host')
    parser.add_argument('-c', '--config-path', required=True, help='specifies projection configuration YAML file path')
    args = parser.parse_args()

    main(args.config_path, args.mount_point, args.data_directory)
    mock_resource.deactivate()
