#!/usr/bin/env python3

import logging
import logging.config
import io
import json
import os
import sys
import time
import xmltodict
import subprocess
from Bio import Entrez


from projections import Projection,  ProjectionDriver, ProjectionTree, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE
from tests.genbank_mock import GenbankMock

logger = logging.getLogger('genbank_projection')


class GenbankDriver(ProjectionDriver):
    def __init__(self, email):
        Entrez.email = email
        Entrez.tool = 'genbank_projection_manager'
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
            esearch_handle = Entrez.esearch(db='nucleotide', term=query[1], retmax=query[2])
            return Entrez.read(esearch_handle)
        else:
            return query[0]

    def get_uri_contents_as_stream(self, query):
        logger.debug('Loading query: %s', query)
        query = query.split(':')
        query_type = query[0]
        if query_type == 'query':
            return Entrez.esearch(db='nuccore', term=query[1], retmax=query[2])
        elif query_type == 'get_gb':
            gb = Entrez.efetch(db='nuccore', id=query[1], rettype='gb', retmode='text')
            return gb.read().encode()
        elif query_type == 'get_fasta':
            fasta = Entrez.efetch(db='nuccore', id=query[1], rettype='fasta', retmode='text')
            with open('efetch_gb.xml', 'w') as f:
                f.write(fasta.read())
            return fasta.read().encode()


class GenbankProjector(Projector):
    def __init__(self, driver, root_projection, root_prototype):
        """
        Initializes SRA Projector with driver, assigns root projection, builds prototype and projection tree.
        :param driver: instance of GenbankDriver
        :param prototype_tree: tree of ProjectionPrototype objects to build projection upon
        """
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree()
        self.root_projection = root_projection
        self.projection_tree.add_projection(self.root_projection, None)

        self.create_projection_tree({'/': root_prototype}, projection_tree=self.projection_tree, parent_projection=self.root_projection)
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

        content = self.driver.get_uri_contents_as_stream(uri)
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

    mock_resource = GenbankMock('http://eutils.ncbi.nlm.nih.gov', 'tests/mock_resource')

    projection_configuration = PrototypeDeserializer('genbank_config.yaml')

    genbank_driver = GenbankDriver('vsvekolkin@parseq.pro')

    root_projection = Projection('/', projection_configuration.root_projection_uri)

    projection_filesystem.projection_manager = GenbankProjector(genbank_driver, root_projection,
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
