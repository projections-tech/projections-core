#!/usr/bin/env python3

import logging
import logging.config
import io
import json
import os
import re
import sys
import time
import xmltodict
import subprocess
from Bio import Entrez
from tests.sra_mock import SRAMock


from projections import Projection,  ProjectionDriver, ProjectionTree, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE


logger = logging.getLogger('sra_projection')


class SRADriver(ProjectionDriver):
    def __init__(self, email):
        """
        Initialize driver telling NCBI using user email and name of program
        :param email:
        :return:
        """
        Entrez.email = email
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

    def load_uri_contents_as_stream(self, uri):
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


class SRAProjector(Projector):
    def __init__(self, driver, root_projection, root_prototype):
        """
        Initializes SRA Projector with driver, assigns root projection, builds prototype and projection tree.
        :param driver: instance of SRADriver
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

        content = self.driver.load_uri_contents_as_stream(uri)
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

    root_projection = Projection('/', projection_configuration.root_projection_uri)

    projection_filesystem.projection_manager = SRAProjector(sra_driver, root_projection,
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
