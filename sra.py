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
from tests.mock import SRAMock


from projections import Projection,  ProjectionDriver, ProjectionTree, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE


logger = logging.getLogger('sra_projection')


class SRADriver(ProjectionDriver):
    def __init__(self, email):
        Entrez.email = email
        Entrez.tool = 'sra_projection_manager'
        self.query_cache = {}

    def get_content(self, query):
        """
        Loads content from SRA using Biopython in driver cache, with queries in format: "query_type:query"
        :param query: str containing query to SRA
        :return: dict of query contents
        """
        query = query.split(':')
        logger.debug('Current query: %s', query)
        query_type = query[0]
        # Query looks as: 'query:Test_species'
        if query_type == 'query':
            # Returns esearch response dict.
            esearch_handle = Entrez.esearch(db='sra', term=query[1], retmax=[2])
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
                self.query_cache[search_query_contents['RUN_SET']['RUN']['@accession']] = search_query_contents['RUN_SET']['RUN']['@accession']+'.sam'
                return search_query_contents
            else:
                return self.query_cache[query[1]]
        else:
            return self.query_cache[query[1]]

    def load_content(self, query):
        logger.debug('Loading query: %s', query)
        query = query.split(':')
        query_type = query[0]

        if query_type == 'get_experiment_runs':
            return subprocess.check_output(['./sratoolkit.2.5.2-ubuntu64/bin/sam-dump', query[1]])
        else:
            return json.dumps(self.query_cache[query[1]]).encode()


class SRAProjector(Projector):
    def __init__(self, driver, root_prototype):
        """
        Initializes SRA Projector with driver, assigns root projection, builds prototype and projection tree.
        :param driver: instance of SRADriver
        :param prototype_tree: tree of ProjectionPrototype objects to build projection upon
        """
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree()
        self.root_projection = Projection('/', 'query:Streptococcus:1')
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
    mock_resource = SRAMock('http://eutils.ncbi.nlm.nih.gov', 'tests/mock_resource')

    projection_configuration = PrototypeDeserializer('sra_config.yaml')

    sra_driver = SRADriver('vsvekolkin@parseq.pro')
    projection_filesystem.projection_manager = SRAProjector(sra_driver, projection_configuration.prototype_tree)
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

