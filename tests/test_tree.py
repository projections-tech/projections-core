import logging
import logging.config
from unittest import TestCase, skip
from projections import Tree

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_tree')

class TestTreeStructure(TestCase):
    def setUp(self):
        self.tree = Tree()
        parent = self.tree
        for i in range(10):
            node = Tree(name=i)
            parent.children[node.name] = node
            parent = node

    def test_find_method(self):
        for i in range(10):
            self.assertTrue(self.tree.find(i, field='name'),
                            msg='Testing Tree find() method for object with name: %s'.format(i))

