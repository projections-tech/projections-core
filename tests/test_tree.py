import logging
import logging.config
from unittest import TestCase, skip
from projections import Tree

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_tree')

class TestTreeStructure(TestCase):
    def setUp(self):
        self.tree = Tree(name='root')
        for i in range(2):
            node = Tree(name=i)
            node.parent = self.tree
            self.tree.children[node.name] = node

        sub_node = node
        for i in range(2,6):
            node = Tree(name=i)
            node.parent = sub_node
            sub_node.children[node.name] = node

        self.last_node = node

    def test_find_method(self):
        """
        Tests Tree find method
        """
        for i in range(6):
            self.assertTrue(self.tree.find(i, field='name'),
                            msg='Testing Tree find() method for object with name: {0}'.format(i))

    def test_path_to_node_method(self):
        """
        Tests Tree path_to_node_method
        """
        # Checking path to root node which is empty
        self.assertListEqual([], self.tree.path_to_node())

        # Checking path to children of root node
        for i in range(2):
            current_node = self.tree.find(i, field='name')
            self.assertListEqual(['root'],
                                 [n.name for n in current_node.path_to_node()],
                                 msg='Checking path to node: {}'.format(i))
        # Checking path to children for subnode with name '1'
        for i in range(2, 6):
            current_node = self.tree.find(i, field='name')
            self.assertListEqual(['root', 1],
                                 [n.name for n in current_node.path_to_node()],
                                 msg='Checking path to node: {}'.format(i))

    def test_node_descendants(self):
        """
        Checks node descendants method of Tree
        """
        expected_child_lists = [[0, 1], [2, 3, 4, 5]]
        for children_dict, exp_list in zip((self.tree.children, self.tree.find(1, field='name').children), expected_child_lists):
            child_list = sorted(children_dict.keys())
            self.assertListEqual(child_list, exp_list, msg='Checking node descendants method.')

    def test_traverse_preorder_method(self):
        """
        Checks order of preorder traversal of Tree
        """
        self.assertListEqual(['root']+list(range(6)),
                             [n.name for n in self.tree.traverse_pre_order()],
                             msg='Checking traverse order for preorder method.')

    def test_traverse_post_order(self):
        """
        Checks order of post order traversal of Tree
        """
        traverse_order = [0, 2, 3, 4, 5, 1, 'root']
        traverse_result = [n.name for n in self.tree.traverse_post_order()]
        self.assertListEqual(traverse_order, traverse_result,
                             msg='Checking order of traverse_post_order method: {}'.format(traverse_result))

    def test_traverse_breadth_first(self):
        """
        Checks order of breadth first traversal of Tree
        """
        traverse_order = ['root', 1, 5, 4, 3, 2, 0]
        traverse_result = [n.name for n in self.tree.traverse_breadth_first()]
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