#!/usr/bin/env python3

import logging
import logging.config
import io
import json
import os
import stat
import sys
import time
import xmltodict
import subprocess
from Bio import Entrez

from projections import Projection,  ProjectionDriver, ProjectionTree, Projector, ProjectionPrototype
from filesystem import ProjectionFilesystem
from fuse import FUSE


logger = logging.getLogger('sra_projection')

class SRADriver(ProjectionDriver):
    def __init__(self, email):
        Entrez.email = email
        Entrez.tool = 'sra_projection_manager'
        self.query_cache = {}

    def get_content(self, query):
        query = query.split(':')
        logger.debug('Current query: %s', query)
        query_type = query[0]

        if query_type == 'query':
            esearch_handle = Entrez.esearch(db='sra', term=query[1], retmax=[2])
            return Entrez.read(esearch_handle)

        elif query_type == 'search_id':
            if query[1] not in self.query_cache:
                fetch_handler = Entrez.efetch(db='sra', id=query[1])
                sample_dict = xmltodict.parse(fetch_handler.read())
                search_query_contents = sample_dict['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
                self.query_cache[query[1]] = search_query_contents
                self.query_cache[search_query_contents['EXPERIMENT']['@accession']] = search_query_contents
                sample_run_set = search_query_contents['RUN_SET']['RUN']
                # Dealing with minor inconsistency here, sometimes RUN field contains run dict, sometimes list of run dicts
                if not type(sample_run_set) == 'list':
                    # Deal with inconsistency by creating list anyway
                    sample_run_set = [sample_run_set]
                for run in sample_run_set:
                    sample_id = run['@accession']
                    self.query_cache[sample_id+'.sam'] = sample_id

                return search_query_contents
            else:
                return self.query_cache[query[1]]
        else:
            return self.query_cache[query[1]]

    def load_content(self, uri):
        pass


class SRAProjectionManager:
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

            projections.append(experiment_projection)
            projections.append(experiment_metadata_projection)

            sample_run_set = experiment_metadata['RUN_SET']['RUN']
            # Dealing with minor inconsistency here, sometimes RUN field contains run dict, sometimes list of run dicts
            if not type(sample_run_set) == 'list':
                # Deal with inconsistency by creating list anyway
                sample_run_set = [sample_run_set]

            for run in sample_run_set:
                sample_id = run['@accession']

                sample_path = os.path.join(experiment_path, sample_id + '.sam')
                sample_projection = Projection(sample_path, 'test')
                sample_projection.sample_id = sample_id
                projections.append(sample_projection)

                logger.debug('Run accession: %s', sample_id)

        for p in projections:
            self.projections[p.path] = p

    def open_resource(self, path):
        projection_on_path = self.projections[path]

        _, resource_file_extension = os.path.splitext(path)
        logger.debug('Resource file extension: %s', resource_file_extension)

        if resource_file_extension == '.json':
            content = projection_on_path.metadata.encode()
        elif resource_file_extension == '.sam':
            content = subprocess.check_output(['./sratoolkit.2.5.2-ubuntu64/bin/sam-dump', projection_on_path.sample_id])

        logger.info('Got path content: %s\n', path)

        self.projections[path].size = len(content)

        file_header = 3
        resource_io = io.BytesIO(content)

        return file_header, resource_io

class SRAProjector(Projector):
    def __init__(self, driver):
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree()
        self.root_projection = Projection('/', 'query:Streptococcus:2')
        self.projection_tree.add_projection(self.root_projection, None)

        prototypes = self.prepare_prototypes()
        self.create_projection_tree(prototypes, projection_tree=self.projection_tree, parent_projection=self.root_projection)
        self.projections = self.projection_tree.projections

    def prepare_prototypes(self):
        query_prototype = ProjectionPrototype('directory')
        query_prototype.name = "environment['QueryTranslation']"
        query_prototype.uri = "['search_id:'+id for id in environment['IdList']]"

        experiment_prototype = ProjectionPrototype('directory', query_prototype)
        experiment_prototype.name = "environment['EXPERIMENT']['@accession']"
        experiment_prototype.uri = "['experiment_id:'+environment['EXPERIMENT']['@accession']]"

        run_prototype = ProjectionPrototype('directory', experiment_prototype)
        run_prototype.name = "'test'"
        run_prototype.uri = "['sample_id:'+environment]"

        query_prototype.children[experiment_prototype.name] = experiment_prototype
        experiment_prototype.children[run_prototype.name] = run_prototype

        return {'/': query_prototype}

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
    sra_driver = SRADriver('vsvekolkin@parseq.pro')
    projection_filesystem.projection_manager = SRAProjector(sra_driver)
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

