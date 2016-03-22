__author__ = 'abragin'

import getpass
import json
import logging
import logging.config
import re
from unittest import TestCase

import psycopg2

from db_projector import DBProjector
from projections import ProjectionPrototype, ProjectionDriver, PrototypeDeserializer

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_projections')


class TestDriver(ProjectionDriver):
    """
    This class does not perform actual testing. It's providing test data and may be replaced with mock object.
    """

    def get_uri_contents_as_dict(self, uri):
        logger.info('Requesting content for uri: %s', uri)

        if not uri:
            logger.info('Returning empty array')
            return []

        if uri == 'experiments':
            with open('tests/json/experiments.json') as f:
                content = json.load(f)
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
            return {'meta': "This is BAM file"}

        assert False is True, 'Test driver can\'t handle resource request, aborting!'


class TestPrototypeDeserializer(TestCase):
    def setUp(self):
        self.deserializer = PrototypeDeserializer('tests/test_projection_config.yaml')

    def test_prototype_deserialization(self):
        """
        Tests test_projection_config deserialization in to ProjectionPrototype tree
        """
        root_prototype = self.deserializer.prototype_tree
        # Test if nodes are ProjectionPrototype instances
        test_nodes_list = [n for n in root_prototype.get_tree_nodes()]
        for n in test_nodes_list:
            self.assertIsInstance(n, ProjectionPrototype,
                                  msg='Checking if object: {0} is instance of ProjectionPrototype'.format(
                                      root_prototype))
        # Test correctness of "name" fields of nodes
        expected_names = ['root_dir', 'results_dir', 'test_bam.bam', 'test_vcf.vcf']
        test_names = [n.name for n in root_prototype.get_tree_nodes()]
        for element in expected_names:
            self.assertIn(element, test_names,
                          msg='Checking existance of projection with name {} in a tree.'.format(element))
        # Test correctness of "uri" fields of nodes
        expected_uri = ['[object["uri"] for object in environment]', "environment['results']",
                        "[environment['data_vcf']]", "[environment['data_bam']]"]
        test_pre_order_uri = [n.uri for n in root_prototype.get_tree_nodes()]
        for element in expected_uri:
            self.assertIn(element, test_pre_order_uri,
                          msg='Checking existance of projection with uri {} in a tree.'.format(element))
        # Test correctness of "type" fields of nodes
        expected_types = ['directory', 'directory', 'file', 'file']
        test_pre_order_uri = [n.type for n in root_prototype.get_tree_nodes()]
        self.assertListEqual(expected_types, test_pre_order_uri, msg='Checking if prototypes types are correct.')


class TestProjector(TestCase):
    """
    Tests DBProjector projection tree building
    """

    @classmethod
    def setUpClass(cls):
        # Initializing database connection which will be used during tests
        cls.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))
        # Creating cursor, which will be used to interact with database
        cls.cursor = cls.db_connection.cursor()

        cls.projection_driver = TestDriver()

        cls.experiment_prototype = ProjectionPrototype('directory')
        cls.experiment_prototype.name = " replace($.displayName, ' ', '_') "
        cls.experiment_prototype.uri = ' $.objects.uri '

        result_prototype = ProjectionPrototype('directory')
        result_prototype.name = " split($.filesystempath, '/')[-1] "
        result_prototype.uri = " $.results "

        bam_prototype = ProjectionPrototype('file')
        bam_prototype.name = " split($.environment.data, '/')[-1] "
        bam_prototype.uri = " $.data "

        result_prototype.parent = cls.experiment_prototype
        cls.experiment_prototype.children[result_prototype.name] = result_prototype
        bam_prototype.parent = result_prototype
        result_prototype.children[bam_prototype.name] = bam_prototype

    @classmethod
    def tearDownClass(cls):
        # Removing test projection entries from projections db
        cls.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_projection' ")
        cls.db_connection.commit()
        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        self.projector = DBProjector('test_projection',
                                     self.projection_driver,
                                     self.experiment_prototype,
                                     'experiments')

    def tearDown(self):
        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM tree_table WHERE projection_name='test_projection' ")
        self.db_connection.commit()

    def test_create_projection_tree(self):
        """
        Testing projection tree creation with projection prototypes.
        """

        dir_paths = ['/', '/experiment_0', '/experiment_1', '/experiment_2',
                     '/experiment_1/result_1', '/experiment_1/result_2',
                     '/experiment_2/result_3', '/experiment_2/result_4', '/experiment_2/result_5']

        created_projections = self.projector.list_projections()

        logger.info('Created projections: %s', created_projections)
        for dir_path in dir_paths:
            logger.info('Checking projection on path: %s', dir_path)

            self.assertIn(dir_path, created_projections, 'Check that projection exists')

            projection_stats = self.projector.get_attributes(dir_path)

            self.assertTrue(projection_stats['st_mode'] == 16895, 'Check that this is a directory projection')

        file_paths = ['/experiment_1/result_1/1.bam', '/experiment_1/result_2/2.bam',
                      '/experiment_2/result_3/3.bam', '/experiment_2/result_4/4.bam', '/experiment_2/result_5/5.bam']

        for file_path in file_paths:
            logger.info('Checking file projection on path: %s', file_path)
            self.assertIn(file_path, created_projections, 'Check that projection exists')

            projection_stats = self.projector.get_attributes(file_path)

            self.assertTrue(projection_stats['st_mode'] == 33279, 'Check that this is a file projection')

    def test_projection_deletion(self):
        """
        Tests projection deletion correctness
        """
        created_projections = self.projector.list_projections()

        # Removing entire directory
        self.projector.remove_projection('/experiment_1')

        current_projections = self.projector.list_projections()

        removed_projections = ['/experiment_1', '/experiment_1/result_1',
                               '/experiment_1/result_1/1.bam', '/experiment_1/result_2',
                               '/experiment_1/result_2/2.bam']

        # Checking if projections where removed properly
        for removed_projection in removed_projections:
            self.assertNotIn(removed_projection, current_projections, 'Checking directory deletion')

        # Testing leaf projections removal
        projections_to_remove = ['/experiment_2/result_4/4.bam', '/experiment_2/result_5/5.bam']
        # Removing leaf projections
        for path in projections_to_remove:
            self.projector.remove_projection(path)

        current_projections = self.projector.list_projections()

        for path in projections_to_remove:
            self.assertNotIn(path, current_projections, 'Checking file deletion')

    def test_projections_move(self):
        """
        Tests projections move correctness
        """
        created_projections = self.projector.list_projections()

        # Moving dir experiment_1 into dir experiment_2
        self.projector.move_projection('/experiment_1', '/experiment_2')

        current_projections = self.projector.list_projections()

        moved_dir_paths = ['/experiment_2/experiment_1', '/experiment_2/experiment_1/result_2',
                           '/experiment_2/experiment_1/result_1',
                           '/experiment_2/experiment_1/result_2/2.bam',
                           '/experiment_2/experiment_1/result_1/1.bam']
        for path in moved_dir_paths:
            self.assertIn(path, current_projections, 'Checing node move')
