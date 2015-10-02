__author__ = 'abragin'

import json
import logging
import logging.config
import re
import stat
from unittest import TestCase, skip

from projections import ProjectionPrototype, Projector, Projection, ProjectionTree

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_projector')

class TestDriver(object):

    def get_content(self, uri):
        logger.info('Requesting content for uri: %s', uri)

        if not uri:
            logger.info('Returning empty array')
            return []

        if uri == 'experiments':
            with open('tests/json/experiments.json') as f:
                content = json.load(f)['objects']
            logger.info('Returning content for uri: %s, content: %s', uri, content)
            return content

        match = re.match('experiments/(\d+)', uri)
        if match:
            id = match.groups()[0]
            logger.info('Requesting experiment data with id: %s', id)
            with open('tests/json/experiment_{}.json'.format(id)) as f:
                content = json.load(f)
            logger.info('Returning content for uri: %s, content: %s', uri, content)
            return content

        match = re.match('results/(\d+)', uri)
        if match:
            id = match.groups()[0]
            logger.info('Requesting result data with id: %s', id)
            return {'id': id, 'content': 'Result of some experiment',
                    'filesystempath': '/tmp/result_{}'.format(id), 'data': 'data/{}.bam'.format(id)}

        match = re.match('data/(\d+).bam', uri)
        if match:
            id = match.groups()[0]
            logger.info('Requesting result data with id: %s', id)
            return "This is BAM file"

        assert False is True, 'Test driver can\'t handle resource request, aborting!'

class TestProjector(TestCase):

    def test_create_projection_tree(self):
        """
        Testing projection tree creation with projection prototypes.

        :return:
        """

        root = Projection('/', 'experiments')

        experiment_prototype = ProjectionPrototype('directory')
        experiment_prototype.name = "content['displayName'].replace(' ', '_')"
        experiment_prototype.uri = '[object["uri"] for object in environment]'

        result_prototype = ProjectionPrototype('directory')
        result_prototype.name = "path.split(content['filesystempath'])[1]"
        result_prototype.uri = "environment['results']"

        bam_prototype = ProjectionPrototype('file')
        bam_prototype.name = "path.split(environment['data'])[1]"
        bam_prototype.uri = "[environment['data']]"

        experiment_prototype.children = [result_prototype]
        result_prototype.children = [bam_prototype]

        # Create projection tree with the prototypes provided
        projector = Projector(TestDriver())
        projection_tree = ProjectionTree()
        projection_tree.add_projection(root, None)

        projector.create_projection_tree([experiment_prototype], projection_tree=projection_tree, parent_projection=root)

        dir_paths = ['/', '/experiment_0', '/experiment_1', '/experiment_2',
                 '/experiment_1/result_1', '/experiment_1/result_2',
                 '/experiment_2/result_3', '/experiment_2/result_4', '/experiment_2/result_5']

        for dir_path in dir_paths:
            logger.info('Checking projection on path: %s', dir_path)
            self.assertTrue(dir_path in projection_tree.projections, 'Check that projection exists')
            projection = projection_tree.projections[dir_path]
            self.assertTrue(projection.type == stat.S_IFDIR, 'Check that this is a directory projection')

        file_paths = ['/experiment_1/result_1/1.bam', '/experiment_1/result_2/2.bam',
                      '/experiment_2/result_3/3.bam', '/experiment_2/result_4/4.bam', '/experiment_2/result_5/5.bam']

        for file_path in file_paths:
            logger.info('Checking file projection on path: %s', file_path)
            self.assertTrue(file_path in projection_tree.projections, 'Check that projection exists')
            projection = projection_tree.projections[file_path]
            self.assertTrue(projection.type == stat.S_IFREG, 'Check that this is a file projection')
