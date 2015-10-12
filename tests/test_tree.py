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

    def test_traverse_preorder_method(self):
        self.assertListEqual(['root']+list(range(6)),
                             [n.name for n in self.tree.traverse_pre_order()],
                             msg='Checking traverse order for preorder method.')

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

    def test_traverse_post_order(self):
        traverse_order = [0, 2, 3, 4, 5, 1, 'root']
        traverse_result = [n.name for n in self.tree.traverse_post_order()]
        self.assertListEqual(traverse_order, traverse_result)

    def test_traverse_breadth_first(self):
        traverse_order = [0, 2, 3, 4, 5, 1, 'root']
        traverse_result = [n.name for n in self.tree.traverse_breadth_first()]
        self.assertListEqual(traverse_order, traverse_result)