#    Copyright 2016  Anton Bragin, Victor Svekolkin
#
#    This file is part of Projections.
#
#    Projections is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Projections is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Projections.  If not, see <http://www.gnu.org/licenses/>.

import logging
import logging.config
import os
import pprint
import threading

import yaml

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('projections')


class Node(object):
    """
    Generic tree class which holds common methods for "tree" structure, must be subclassed to create Prototype and
    Projection classes
    """

    def __init__(self, name=None, data=None):
        """
        Initializes Node with optional arguments name and data
        :param name: Name of current node
        :param data: Field that holds data associated with current node
        """
        self.name = name
        self.data = data
        self.root = None
        self.parent = None
        self.children = {}
        # Be aware of atomicity of operations
        self.lock = threading.Lock()

    def __split_path(self, path):
        """
        Split path string in parts preserving path root
        :param path: path string
        :return: list of path parts strings
        """
        temp = []
        while True:
            head, tail = os.path.split(path)
            if head == '/' or head == '':
                if tail == '':
                    return [head]
                else:
                    temp.append(tail)
                    if not head == '':
                        temp.append(head)
                    return temp[::-1]
            temp.append(tail)
            path = head

    def add_child(self, node):
        """
        Adds node to current node child nodes dict
        :param node: Node object instance
        """
        self.lock.acquire()
        node.parent = self
        self.children[node.name] = node
        self.lock.release()

    def get_parent_nodes(self):
        """
        Returns nodes list from root to current node
        """
        result = []
        parent, child = self.parent, self
        while parent:
            result.append(parent)
            parent, child = parent.parent, parent
        return result[::-1]

    def get_path(self):
        """
        Returns path to node from root string
        :return: path to node from root string
        """
        parents = self.get_parent_nodes()

        # Node that have no parents is root node, so it`s name is '/'
        if parents:
            return ''.join(['/'] + ['/'.join([n.name for n in parents[1:]] + [self.name])])
        else:
            return '/'

    def get_children(self):
        """
        Returns children iterator of current node
        """
        return self.children.values()

    def get_children_names(self):
        """
        Returns iterator of children names
        :return: iterator of children names strings
        """
        return self.children.keys()

    def find_node_by_path(self, path_to_node):
        """
        Finds node in a tree according to path. Path to node is relative
        """
        path = self.__split_path(path_to_node)
        temp_node = self
        # If node have parent start check of node name from second element of path, first always being "/"
        if temp_node.parent:
            head, tail = path[1], path[2:]
        else:
            head, tail = path[0], path[1:]

        if temp_node.name == head:
            for item in tail:
                if item in temp_node.children:
                    temp_node = temp_node.children[item]
                else:
                    return None
        else:
            return None
        return temp_node

    def get_tree_nodes(self):
        """
        Generator that traverses current tree, yields nodes in current tree, level by level
        """
        yield self
        for c in self.children.values():
            for v in c.get_tree_nodes():
                yield v

    def remove_node_by_path(self, path):
        """
        Removes node on path from current tree
        """
        try:
            node_to_remove = self.find_node_by_path(path)
            if node_to_remove:
                logger.debug('Removing node: %s', node_to_remove.name)
                # Dict pop method is considered atomic
                node_to_remove.parent.children.pop(node_to_remove.name)
        except:
            raise KeyError('Attempting to remove node that is not in a tree!')

    def __str__(self):
        return pprint.pformat([n.get_path() for n in self.get_tree_nodes()])


class ProjectionPrototype(Node):
    """
    The class objects describe nodes in projection logical structure.
    Every Prototype object may have 0..many projections associated with it.

    Serialized form of the ProjectionPrototype objects hierarchy is THE way to describe projections.
    """

    def __init__(self, type, parent=None, name=None):
        """
        Create ProjectionPrototype object.

        :param type: describe the type of the generated projections. Current implementation uses 'directory' and 'file' types
        """
        # Initialize Node class, passing current object in Node data field
        super().__init__(name, self)
        self.parent = parent
        self.type = type
        self.meta_link = [None]
        # TODO: consider logical synchronization of name and uri
        # Dialect specific description that is used as a generator for projection uri's
        self.uri = None
        self.context = None

    def get_context(self):
        """
        Get context of current node defined by contexts of upper nodes
        """
        return [node.context for node in self.get_parent_nodes()]

    def __str__(self):
        return 'ProjectionPrototype for {}'.format(self.uri)

    def __repr__(self):
        return self.__str__()


class PrototypeDeserializer(object):
    """
    Used to deserialize Prototype tree from YAML configuration files
    """

    def __init__(self, data_path):
        """
        Initialize class, passing path to YAML configuration file
        """
        yaml_stream = self.read_yaml_file(data_path)
        self.prototype_tree, self.resource_uri, self.root_projection_uri, self.driver_config_path = self.read_projections(
            yaml_stream)

    def get_prototypes_tree(self, yaml_dict, parent=None):
        """
        Constructs ProjectionPrototype tree from root dictionary
        :param yaml_dict: dict for root element
        :param parent: parent prototype for tree node
        :return root ProjectionPrototype object
        """
        pp = ProjectionPrototype(type=yaml_dict['type'])
        pp.name = yaml_dict['name']
        pp.uri = yaml_dict['uri']
        pp.parent = parent
        if 'meta_link' in yaml_dict:
            pp.meta_link = yaml_dict['meta_link']

        if isinstance(yaml_dict['children'], dict):
            pp.children = {x[0]: self.get_prototypes_tree(x[1], parent=pp) for x in yaml_dict['children'].items()}
            return pp
        else:
            return pp

    def read_yaml_file(self, file_path):
        """
        Opens yaml file on path
        :param file_path: path to yaml file string
        :return: yaml file stream
        """
        return open(file_path)

    def read_projections(self, yaml_stream):
        """
        Read YAML configuration file on data_path and construct Prototype tree from it
        :param yaml_stream yaml file stream
        :return root ProjectionPrototype object
        """
        with yaml_stream:
            yaml_dict = yaml.safe_load(yaml_stream)

        resource_uri = yaml_dict['resource_uri']
        root_projection_uri = yaml_dict['root_projection_uri']
        driver_config_path = yaml_dict['driver_config_path']

        return self.get_prototypes_tree(yaml_dict['root']), resource_uri, root_projection_uri, driver_config_path


class ProjectionDriver(object):
    """
    Object that has get_uri_content(uri) method returning array of Python dictionaries.
    """

    def read_config(self, daemon_script_dir, config_path):
        # Opening driver configuration
        with open(os.path.join(daemon_script_dir, config_path)) as yaml_stream:
            return yaml.safe_load(yaml_stream)

    def get_uri_contents_as_dict(self, uri):
        raise NotImplemented('Implement metadata retrieval from some projection backend.')

    def get_uri_contents_as_bytes(self, uri):
        raise NotImplemented('Implement data stream retrieval from some projection backend.')
