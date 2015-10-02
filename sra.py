#!/usr/bin/env python3

import logging
import io
import json
import os
import stat
import sys
import time
import xmltodict
import urllib.request
from urllib.parse import urljoin
from Bio import Entrez

from projections import Projection, ProjectionManager
from filesystem import ProjectionFilesystem
from fuse import FUSE


logger = logging.getLogger('sra_projection')



class SRAProjectionManager(ProjectionManager):
    def __init__(self, email, query, num_of_results=1):
        self.query = query
        self.num_of_results = num_of_results
        logger.info('Creating SRA projection for query: %s', self.query)
        self.setup_biopython(email)

        # TODO: switch to tree-like structure instead of manual path parsing
        self.projections = {}
        self.create_projections()

    def setup_biopython(self, email):
        Entrez.email = email
        Entrez.tool = 'sra_projection_manager'

    def create_projections(self):
        """
        Creates projections for SRA search request
        :return: list of projections
        """
        projections = []

        self.search_handle = Entrez.esearch(db='sra', term=self.query, retmax=self.num_of_results)
        self.search_results = Entrez.read(self.search_handle)

        query_base_path = '/' + self.query
        query_projection = Projection(query_base_path, 'test')
        query_projection.type = stat.S_IFDIR
        projections.append(query_projection)

        results_ids = self.search_results['IdList']
        for ids in results_ids:
            logger.debug('Query ID: %s', ids)
            fetch_handler = Entrez.efetch(db='sra', id=ids)
            sample_dict = xmltodict.parse(fetch_handler.read())

            experiment_metadata = sample_dict['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
            experiment_id = experiment_metadata['EXPERIMENT']['@accession']

            experiment_path = os.path.join(query_base_path, experiment_id)
            experiment_projection = Projection(experiment_path, '')
            experiment_projection.type = stat.S_IFDIR

            ex_metadata_path = os.path.join(experiment_path, 'experiment_metadata.json')
            ex_metadata_uri = 'www.ncbi.nlm.nih.gov/sra/{0}'.format(experiment_id)
            experiment_metadata_projection = Projection(ex_metadata_path, ex_metadata_uri)
            experiment_metadata_projection.metadata = json.dumps(experiment_metadata)
            logger.debug('Experiment metadata: %s', experiment_metadata_projection.metadata)

            projections.append(experiment_projection)
            projections.append(experiment_metadata_projection)

            sample_run_set = experiment_metadata['RUN_SET']['RUN']
            sample_id = sample_run_set['@accession']

            sample_path = os.path.join(experiment_path, sample_id+'.bam')
            sample_projection = Projection(sample_path, 'test')
            projections.append(sample_projection)

            logger.debug('Run accession: %s', sample_id)

        for p in projections:
            self.projections[p.path] = p

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
        projection_on_path = self.projections[path]

        _, resource_file_extension = os.path.splitext(path)

        logger.debug('Resource file extension: %s', resource_file_extension)

        if resource_file_extension == '.json':
            content = projection_on_path.metadata.encode()
        elif resource_file_extension == '.bam':
            content = b'Test BAM!'

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
    projection_filesystem.projection_manager = SRAProjectionManager('vsvekolkin@parseq.pro', 'Streptococcus', 1)
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

