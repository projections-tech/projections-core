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
        self.tree = Tree(name=None)

        self.first_level_names = []
        self.second_level_names = []
        self.third_level_names = []

        for i in range(3):
            f_name = '{0}'.format(i)
            self.first_level_names.append(f_name)
            first_level = Tree(name=f_name)
            self.tree.add_child(first_level)
            for j in range(3):
                s_name = '{0}.{1}'.format(i, j)
                self.second_level_names.append(s_name)
                second_level = Tree(s_name)
                first_level.add_child(second_level)
                for k in range(3):
                    t_name = '{0}.{1}.{2}'.format(i, j, k)
                    self.third_level_names.append(t_name)
                    third_level = Tree(name=t_name)
                    second_level.add_child(third_level)

    def test_get_path(self):
        """
        Tests Tree get_path
        """
        # Checking get_path for root node which is empty list
        self.assertListEqual([], self.tree.get_path())

        for first_level_child in self.tree.get_children():
            self.assertListEqual([None],
                                 [n.name for n in first_level_child.get_path()],
                                 msg='Checking path to first level nodes of a tree')

            for second_level_child in first_level_child.get_children():
                self.assertListEqual([None, first_level_child.name],
                                 [n.name for n in second_level_child.get_path()],
                                 msg='Checking path to second level nodes of a tree')

                for third_level_child in second_level_child.get_children():
                    self.assertListEqual([None, first_level_child.name, second_level_child.name],
                                     [n.name for n in third_level_child.get_path()],
                                     msg='Checking path to third level nodes of a tree')

    def test_get_children(self):
        """
        Checks get_children method of Tree
        """
        for first_level_child in self.tree.get_children():
            self.assertIn(first_level_child.name, self.first_level_names,
                          msg='Checking name to first level nodes of a tree')

            for j, second_level_child in enumerate(first_level_child.get_children()):
                self.assertIn(second_level_child.name, self.second_level_names,
                                 msg='Checking name to second level nodes of a tree')

                for k, third_level_child in enumerate(second_level_child.get_children()):
                    self.assertIn(third_level_child.name, self.third_level_names,
                                  msg='Checking name to third level nodes of a tree')

    def test_find_node_by_path(self):
        """
        Checks if node is on path using Tree find_node_by_path method
        """
        tree = Tree(name='/')
        node_1 = Tree(name='experiments')
        node_1_1 = Tree(name='a')
        node_1_2 = Tree(name='b')
        node_2 = Tree(name='results')
        node_2_1 = Tree(name='a')
        node_2_2 = Tree(name='b')

        tree.add_child(node_1)
        node_1.add_child(node_2)
        node_1.add_child(node_1_1)
        node_1.add_child(node_1_2)
        node_2.add_child(node_2_1)
        node_2.add_child(node_2_2)

        path_to_node_name_dict = {'/': '/', '/experiments': 'experiments', '/experiments/results': 'results',
                                  '/experiments/a': 'a', '/experiments/b': 'b', '/experiments/results/a': 'a',
                                  '/experiments/results/b': 'b'}

        for path, node_name in path_to_node_name_dict.items():
            node_by_path = tree.find_node_by_path(path)
            self.assertTrue(node_by_path.name == node_name,
                            msg='Checking if node is on path: {0}'.format(path))

        # Test find node by path starting not from root node
        path = 'experiments/results'
        node_by_path = node_1.find_node_by_path(path)
        self.assertTrue(node_by_path.name == 'results',
                        msg='Checking if node: {0} is on path: {1}'.format(node_by_path.name, path))


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
                                  msg='Checking if object: {0} is instance of ProjectionPrototype'.format(root_prototype))
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