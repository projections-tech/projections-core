__author__ = 'abragin'

import json
import logging
import logging.config
import re
import stat
import json

from unittest import TestCase, skip

from projections import Tree, ProjectionPrototype, Projector, Projection, ProjectionTree, ProjectionDriver, PrototypeDeserializer

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_projections')


class TestDriver(ProjectionDriver):
    """
    This class does not perform actual testing. It's providing test data and may be replaces with mock object.
    """
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

        """

        # Root projection has ho associated prototypes.
        # This behavior may be changed to provide more uniform approach.
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

        result_prototype.parent = experiment_prototype
        experiment_prototype.children[result_prototype.name] = result_prototype
        bam_prototype.parent = result_prototype
        result_prototype.children[bam_prototype.name] = bam_prototype

        # Create projection tree with the prototypes provided
        projector = Projector(TestDriver())
        projection_tree = ProjectionTree()
        projection_tree.add_projection(root, None)

        projector.create_projection_tree({'/': experiment_prototype}, projection_tree=projection_tree, parent_projection=root)

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


class TestTree(TestCase):
    def setUp(self):
        self.tree = Tree(name='root')
        for i in range(3):
            first_level = Tree(name='{0}'.format(i))
            first_level.parent = self.tree
            self.tree.add_child(first_level)
            for j in range(3):
                second_level = Tree(name='{0}.{1}'.format(i, j))
                second_level.parent = first_level
                first_level.add_child(second_level)
                for k in range(3):
                    third_level = Tree(name='{0}.{1}.{2}'.format(i, j, k))
                    third_level.parent = second_level
                    second_level.add_child(third_level)

    def test_find(self):
        """
        Tests Tree find method
        """
        for i in range(3):
            self.assertIsInstance(self.tree.find('{0}'.format(i)), Tree,
                                  msg='Testing Tree find() method for object with name: {0}'.format(i))
            for j in range(3):
                self.assertIsInstance(self.tree.find('{0}.{1}'.format(i, j)), Tree,
                                      msg='Testing Tree find() method for object with name: {0}.{1}'.format(i, j))
                for k in range(3):
                    self.assertIsInstance(self.tree.find('{0}.{1}.{2}'.format(i, j, k)), Tree,
                                          msg='Testing Tree find() method for object with name: {0}.{1}.{2}'.format(i, j, k))

    def test_get_path(self):
        """
        Tests Tree get_path
        """
        # Checking get_path for root node which is empty list
        self.assertListEqual([], self.tree.get_path())

        for i, first_level_child in enumerate(self.tree.get_children()):
            self.assertListEqual(['root'],
                                 [n.name for n in first_level_child.get_path()],
                                 msg='Checking path to first level nodes of a tree')

            for j, second_level_child in enumerate(first_level_child.get_children()):
                self.assertEqual(['root', str(i)],
                                 [n.name for n in second_level_child.get_path()],
                                 msg='Checking path to second level nodes of a tree')

                for k, third_level_child in enumerate(second_level_child.get_children()):
                    self.assertEqual(['root', str(i),'{0}.{1}'.format(i, j)],
                                     [n.name for n in third_level_child.get_path()],
                                     msg='Checking path to third level nodes of a tree')

    def test_get_children(self):
        """
        Checks get_children method of Tree
        """
        for i, first_level_child in enumerate(self.tree.get_children()):
            self.assertEqual(str(i),
                             first_level_child.name,
                             msg='Checking name to first level nodes of a tree')

            for j, second_level_child in enumerate(first_level_child.get_children()):
                self.assertEqual('{0}.{1}'.format(i,j),
                                 second_level_child.name,
                                 msg='Checking name to second level nodes of a tree')

                for k, third_level_child in enumerate(second_level_child.get_children()):
                    self.assertEqual('{0}.{1}.{2}'.format(i, j, k),
                                     third_level_child.name,
                                     msg='Checking name to third level nodes of a tree')

    def test_traverse_preorder(self):
        """
        Checks order of preorder traversal of Tree
        """
        traverse_order = ['root', '0', '0.0', '0.0.0', '0.0.1', '0.0.2', '0.1', '0.1.0', '0.1.1', '0.1.2', '0.2',
                          '0.2.0', '0.2.1', '0.2.2', '1', '1.0', '1.0.0', '1.0.1', '1.0.2', '1.1', '1.1.0', '1.1.1',
                          '1.1.2', '1.2', '1.2.0', '1.2.1', '1.2.2', '2', '2.0', '2.0.0', '2.0.1', '2.0.2', '2.1',
                          '2.1.0', '2.1.1', '2.1.2', '2.2', '2.2.0', '2.2.1', '2.2.2']
        self.assertListEqual(traverse_order,
                             [n.name for n in self.tree.traverse_pre_order()],
                             msg='Checking traverse order for preorder method.')

    def test_traverse_post_order(self):
        """
        Checks order of post order traversal of Tree
        """
        traverse_order = ['0.0.0', '0.0.1', '0.0.2', '0.0', '0.1.0', '0.1.1', '0.1.2', '0.1', '0.2.0', '0.2.1', '0.2.2',
                          '0.2', '0', '1.0.0', '1.0.1', '1.0.2', '1.0', '1.1.0', '1.1.1', '1.1.2', '1.1', '1.2.0', '1.2.1',
                          '1.2.2', '1.2', '1', '2.0.0', '2.0.1', '2.0.2', '2.0', '2.1.0', '2.1.1', '2.1.2', '2.1', '2.2.0',
                          '2.2.1', '2.2.2', '2.2', '2', 'root']
        traverse_result = [n.name for n in self.tree.traverse_post_order()]
        self.assertListEqual(traverse_order, traverse_result,
                             msg='Checking order of traverse_post_order method: {}'.format(traverse_result))

    def test_traverse_depth_first(self):
        """
        Checks order of depth first traversal of Tree
        """
        traverse_order = ['root', '2', '2.2', '2.2.2', '2.2.1', '2.2.0', '2.1', '2.1.2', '2.1.1', '2.1.0', '2.0',
                          '2.0.2', '2.0.1', '2.0.0', '1', '1.2', '1.2.2', '1.2.1', '1.2.0', '1.1', '1.1.2', '1.1.1',
                          '1.1.0', '1.0', '1.0.2', '1.0.1', '1.0.0', '0', '0.2', '0.2.2', '0.2.1', '0.2.0', '0.1',
                          '0.1.2', '0.1.1', '0.1.0', '0.0', '0.0.2', '0.0.1', '0.0.0']
        traverse_result = [n.name for n in self.tree.traverse_depth_first()]
        self.assertListEqual(traverse_order, traverse_result,
                             msg='Checking order of traverse_breadth_first: {}'.format(traverse_result))

    def test_find_node_by_path(self):
        """
        Checks if node is on path using Tree find_node_by_path method
        """
        tree = Tree(name='/')
        node_1 = Tree(name='experiments')
        node_1.parent = tree
        node_2 = Tree(name='results')
        node_2.parent = node_1

        tree.children[node_1.name] = node_1
        node_1.children[node_2.name] = node_2

        path = '/experiments/results'
        node_by_path = tree.find_node_by_path(path)
        self.assertTrue(node_by_path.name == 'results',
                        msg='Checking if node: {0} is on path: {1}'.format(node_by_path.name, path))



class TestPrototypeDeserializer(TestCase):
    def setUp(self):
        self.deserializer = PrototypeDeserializer('tests/test_projection_config.yaml')

    def test_prototype_deserialization(self):
        """
        Tests test_projection_config deserialization in to ProjectionPrototype tree
        :return:
        """
        root_prototype = self.deserializer.prototype_tree
        # Test if nodes are ProjectionPrototype instances
        test_pre_order = [n for n in root_prototype.traverse_pre_order()]
        for n in test_pre_order:
            self.assertIsInstance(n, ProjectionPrototype,
                                  msg='Checking if object: {0} is instance of ProjectionPrototype'.format(root_prototype))
        # Test correctness of "name" fields of nodes
        expected_names = ['root_dir', 'results_dir', 'test_bam.bam', 'test_vcf.vcf']
        test_pre_order_names = [n.name for n in root_prototype.traverse_pre_order()]
        for element in expected_names:
            self.assertIn(element, test_pre_order_names,
                          msg='Checking existance of projection with name {} in a tree.'.format(element))
        # Test correctness of "uri" fields of nodes
        expected_uri = ['[object["uri"] for object in environment]', "environment['results']",
                        "[environment['data_vcf']]", "[environment['data_bam']]"]
        test_pre_order_uri = [n.uri for n in root_prototype.traverse_pre_order()]
        for element in expected_uri:
            self.assertIn(element, test_pre_order_uri,
                          msg='Checking existance of projection with uri {} in a tree.'.format(element))
        # Test correctness of "type" fields of nodes
        expected_types = ['directory', 'directory', 'file', 'file']
        test_pre_order_uri = [n.type for n in root_prototype.traverse_pre_order()]
        self.assertListEqual(expected_types, test_pre_order_uri, msg='Checking if prototypes types are correct.')