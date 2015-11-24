__author__ = 'abragin'

import io
import json
import logging
import os
import stat
import time
import threading
import yaml
import pprint
import objectpath
import types
import copy
import itertools

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


class Projection(object):
    """
    Class for holding projection data.

    This is designed to be implementation independent.
    The objects of the class should be logically immutable.

    The purpose of class objects is to store data related to filesystem object representation (name, type and size)
        and mapping of path to resource uri.
    """

    def __init__(self, name, uri, type=stat.S_IFDIR, size=4096):
        """
        Create projection object.

        :param name: name of projection
        :param uri: resource identifier projection is pointing to
        :param type: type of filesystem object. Current implementation supports S_IFDIR and S_IFREG only.
        :param size: size of filesystem object in bytes.
        """
        self.name = name
        self.uri = uri
        self.type = type
        # Should be nonzero to trigger projection object reading
        self.size = size
        self.metadata_uri = None

    def __str__(self):
        return 'Projection from {} to {}'.format(self.uri, self.name)


class ProjectionTree(Node):
    """
    Class for holding tree-like projection structure.
    """

    def __init__(self, p_name, p_uri, p_parent=None, p_type=stat.S_IFDIR, p_size=4096, driver=None):
        super().__init__(name=p_name, data=Projection(name=p_name, uri=p_uri, type=p_type, size=p_size))
        self.parent = p_parent
        self.driver = driver

    def get_projection(self, path):
        """
        Get projection for a path provided.
        :param path: path the projection resides on
        :return: projection or None
        """
        node_on_path = self.find_node_by_path(path)
        return node_on_path.data if node_on_path else None

    def is_managing_path(self, path):
        if self.get_projection(path):
            return True
        else:
            return False

    def get_projections_on_path(self, path):
        """
        Get list of projections on path
        :param path: path string
        :return: list of Projection objects
        """
        logger.info('Requesting projections for path: %s', path)
        node = self.find_node_by_path(path)
        projections = [n.data.name for n in node.get_children()]
        logger.info('Returning projections: %s', projections)
        return projections

    def get_attributes(self, path):
        """
        Get attributes of projection on given path
        """

        assert self.get_projection(path) is not None

        projection = self.get_projection(path)

        now = time.time()
        attributes = dict()

        # Set projection attributes

        # This is implementation specific and should be binded to projector data
        attributes['st_atime'] = now
        # This may be implemented as last projection cashing time is casing is enabled
        attributes['st_mtime'] = now
        # On Unix this is time for metadata modification we can use the same conception
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
        """
        Opens resource on path and returns it`s header and contents stream
        :param path path string
        :return file_header
        :return resource_io
        """
        projection_on_path = self.get_projection(path)
        uri = projection_on_path.uri

        content = self.driver.get_uri_contents_as_bytes(uri)
        logger.info('Got path content: %s\n', path)

        projection_on_path.size = len(content)

        file_header = 3
        resource_io = io.BytesIO(content)

        return file_header, resource_io


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
        self.prototype_tree, self.resource_uri, self.root_projection_uri = self.read_projections(yaml_stream)

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
        return self.get_prototypes_tree(yaml_dict['root']), yaml_dict['resource_uri'], yaml_dict['root_projection_uri']


class ProjectionDriver(object):
    """
    Object that has get_uri_content(uri) method returning array of Python dictionaries.
    """

    def get_uri_contents_as_dict(self, uri):
        raise NotImplemented('Implement metadata retrieval from some projection backend.')

    def get_uri_contents_as_bytes(self, uri):
        raise NotImplemented('Implement data stream retrieval from some projection backend.')


class Projector:
    """
    Class that creates Projection objects and assembles them in ProjectionTree object using Prototypes object.
    """

    def __init__(self, driver, root_projection_uri, prototype_tree):
        """
        :param driver: instance of ProjectionDriver
        :param prototype_tree: tree of ProjectionPrototype objects to build projection upon
        """
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver
        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree(p_name='/', p_uri=root_projection_uri, driver=driver)
        self.create_projection_tree({'/': prototype_tree}, self.projection_tree)

    def fetch_context(self, uri):
        """
        Used to fetch contents of uri in microcode
        :param uri: URI string
        :return: uri contents
        """
        return self.driver.get_uri_contents_as_dict(uri)

    def create_projection_tree(self, prototypes, projection_tree):
        """
        Creates projection tree for a given collection of prototypes.

        Method signature is influenced by recursion implementation and is subject to change.

        :param prototypes: collection of prototypes that are attached to parent projection's prototype
        :param projection_tree: projection tree to extend. Should not be None
        :return: void. The provided projection tree object expected to be used.
        """
        logger.info('Creating projection tree with a prototypes: %s starting from: %s',
                    prototypes, projection_tree.data)
        # Dictionaries that act as a evaluation context for projections. environment is available for both directory and
        # file prototypes, while content for directory prototypes only
        environment = None
        content = None

        metadata_uri = projection_tree.data.uri

        path = os.path

        # This is environment in which projections are created (parent_projection content)
        # TODO: in many cases it means double request to parent projection resource so it should be optimized
        # We don`t want to change driver contents, hence we made deep copy of dict
        environment = copy.deepcopy(self.driver.get_uri_contents_as_dict(projection_tree.data.uri))

        logger.info(' '.join([p.type for p in prototypes.values()]))

        logger.info('Starting prototype creation in the context of resource with uri: %s', projection_tree.data.uri)

        # For every prototype in collection try to create corresponding projections
        for key, prototype in prototypes.items():
            # Set current prototype context to current environment for children node to use
            prototype.context = environment
            # Get context of current node from contexts of parent nodes
            context = prototype.get_context()
            context = context[::-1]

            logger.info('Creating projections for a prototype: %s', prototype)
            # Adding context of upper level prototypes for lower level projections to use
            environment['context'] = context

            # Creating tree of environment contents which will be parsed by ObjectPath
            tree = objectpath.Tree(environment)
            URIs = tree.execute(prototype.uri)

            # Object path sometimes returns generator if user uses selectors, for consistency expand it using
            # list comprehension
            if isinstance(URIs, types.GeneratorType):
                URIs = [el for el in URIs]
            # Treating URIs as list for consistency
            if not isinstance(URIs, list):
                URIs = [URIs]
            logger.info('Prototype %s has projections on URIs: %s', prototype, URIs)

            # We get projection URIs based on environment and prototype properties
            # Every URI corresponds to projection object
            for uri in URIs:
                if prototype.type == 'metadata':
                    # Metadata must be valid JSON file, in other case do conversion
                    environment['file_metadata'] = json.loads(self.driver.get_uri_contents_as_bytes(uri).decode())

                # Get content for a projection
                # We don`t want to change driver contents, hence we made deep copy
                content = copy.deepcopy(self.driver.get_uri_contents_as_dict(uri))
                logger.debug('ENV: %s, CONTENT: %s', environment, content)

                # Adding environment to use by prototype
                content['environment'] = environment
                content['context'] = context

                # Creating tree which will be parsed by ObjectPath
                tree = objectpath.Tree(content)
                name = tree.execute(prototype.name)

                if prototype.type == 'metadata':
                    if isinstance(name, types.GeneratorType):
                        name = [el for el in name]
                    logger.info('NAme contents: %s', name)
                    if not name:
                        break

                # Object path sometimes returns generator if user uses selectors, for consistency expand it using
                # list comprehension
                if isinstance(name, types.GeneratorType):
                    name = [el for el in name]

                child_tree = ProjectionTree(p_name=name, p_uri=uri)
                # Setting projection metadata URI
                child_tree.data.metadata_uri = metadata_uri

                # Add newly created projection to projection tree if prototype is not transparent or metadata
                if prototype.type != 'transparent' and prototype.type != 'metadata':
                    projection_tree.add_child(child_tree)
                    logger.info('Projection created: %s', child_tree)

                if prototype.type == 'directory':
                    # If prototype has children, continue tree building
                    if prototype.children:
                        logger.info('Starting attached prototype projection creation for prototype: %s with children: %s',
                                    prototype, prototype.children)
                        self.create_projection_tree(prototype.children, child_tree)

                elif prototype.type == 'transparent' or prototype.type == 'metadata':
                    # If projection is transparent or metadata, continue prototype tree building passing parent
                    # projection tree and evaluated uri of prototype to it. This behaviour allows to build flat
                    # projections there all prototypes are on one level due to passed common parent projection
                    if prototype.children:
                        projection_tree.data.uri = child_tree.data.uri
                        self.create_projection_tree(prototype.children, projection_tree)

                elif prototype.type == 'file':
                    # NOTE: content variable is not accessible during file projection creation!
                    # This is the point where metadata can be extracted from the file
                    #   but file content should be accessed in this case
                    child_tree.data.type = stat.S_IFREG
                    child_tree.data.size = 1
