__author__ = 'abragin'

import io
import logging
import os
import stat
import time
import threading
import yaml
import pprint

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
        node.parent = self
        self.children[node.name] = node

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
        logger.debug('Splitted path: %s', path)
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

    def __str__(self):
        return pprint.pformat([n.get_path() for n in self.get_tree_nodes()])


class Projection(object):
    """
    Class for holding projection data.

    This is designed to be implementation independent.
    The objects of the class should be logically immutable.

    The purpose of class objects is to store data related to filesystem object representation (name, type and size)
        and mapping of path to resource uri.

    TODO: current implementation is based on path semantics. It should be reconsidered to use name semantics when
    Projection objects holds data on its name only and ProjectionTree object holds path/context data in the form
    of parent/child relationship.

    """

    def __init__(self, path, uri, type=stat.S_IFDIR, size=4096):
        """
        Create projection object.

        :param path: specify path relative to mount point. Path is always prefixed with root '/' sign
        :param uri: resource identifier projection is pointing to
        :param type: type of filesystem object. Current implementation supports S_IFDIR and S_IFREG only.
        :param size: size of filesystem object in bytes.
        :return:
        """
        # TODO: replace 'path' semantics with 'name' semantic
        self.path = path
        self.uri = uri
        self.type = type
        # Should be nonzero to trigger projection object reading
        self.size = size

    def __str__(self):
        return 'Projection from {} to {}'.format(self.uri, self.path)


class ProjectionTree(object):
    """
    Class for holding tree-like projection structure.
    """

    def __init__(self, p_name, p_uri, p_type=stat.S_IFDIR, p_size=4096):
        self.root = Node(name=p_name)
        self.root.data = self
        self.projection = Projection(path=p_name, uri=p_uri, type=p_type, size=p_size)
        # Be aware of atomicity of operations
        self.lock = threading.Lock()

    def get_projection(self, path):
        """
        Get projection for a path provided.
        :param path: path the projection resides on
        :return: projection or None
        """
        node_on_path = self.root.find_node_by_path(path)
        return node_on_path.data if node_on_path else None

    def get_children(self, path):
        """
        Return first-order descendants of the projection on the path specified.
        Implements 'get-directory-content-like' functionality.

        :param path: path to retrieve child projections for
        :return: list of child projections that can be empty
        """
        node_on_path = self.root.find_node_by_path(path)
        if node_on_path:
            return [c.data for c in node_on_path.get_children()]
        else:
            return None

    def add_child(self, proj_tree):
        """
        Adds projection to the current tree.

        :param projection: projection object to add.
        :param parent: parent projection the object should be attached to.
        """
        # This should be done atomically regardless of implementation details.
        self.lock.acquire()

        self.root.add_child(proj_tree.root)

        self.lock.release()

    def __str__(self):
        return 'Projection from {} to {}'.format(self.projection.uri, self.projection.path)


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
        # Dialect specific description that is used as a generator for projection usi's
        self.uri = None
        self.context = None

    def __str__(self):
        return 'ProjectionPrototype for {}'.format(self.uri)

    def __repr__(self):
        return self.__str__()

    def get_context(self):
        """
        Get context of current node defined by contexts of upper nodes
        """
        return [node.context for node in self.get_parent_nodes()]


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

    def get_uri_contents_as_stream(self, uri):
        raise NotImplemented('Implement data stream retrieval from some projection backend.')


class Projector:
    """
    Class that creates Projection object and assembles them in ProjectionTree object using Prototypes object.
    """

    def __init__(self, driver, root_projection_uri, prototype_tree):
        """
        :param driver: instance of ProjectionDriver
        :param prototype_tree: tree of ProjectionPrototype objects to build projection upon
        """
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

        # Initializing projection tree with root projection.
        self.projection_tree = ProjectionTree(p_name='/', p_uri=root_projection_uri)

        self.create_projection_tree({'/': prototype_tree}, projection_tree=self.projection_tree)

    def is_managing_path(self, path):
        if self.projection_tree.get_projection(path):
            return True
        else:
            return False

    def get_projections(self, path):
        logger.info('Requesting projections for path: %s', path)
        projections = [c.projection.path for c in self.projection_tree.get_children(path)]
        logger.info('Returning projections: %s', projections)
        return projections

    def get_attributes(self, path):
        assert self.projection_tree.get_projection(path) is not None

        projection = self.projection_tree.get_projection(path).projection

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
        projection_on_path = self.projection_tree.get_projection(path).projection
        uri = projection_on_path.uri

        content = self.driver.get_uri_contents_as_stream(uri)
        logger.info('Got path content: %s\n', path)

        projection_on_path.size = len(content)

        file_header = 3
        resource_io = io.BytesIO(content)

        return file_header, resource_io

    def fetch_context(self, uri):
        """
        Used to fetch contents of uri in microcode
        :param uri: URI
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
                    prototypes, projection_tree.projection)

        # Dictionaries that act as a evaluation context for projections. environment is available for both directory and
        #   file prototypes, while content for directory prototypes only
        environment = None
        content = None
        fetch_context = self.fetch_context
        path = os.path

        # This is environment in which projections are created (parent_projection content)
        # TODO: in many cases it means double request to parent projection resource so it should be optimized

        environment = self.driver.get_uri_contents_as_dict(projection_tree.projection.uri)
        logger.info('Starting prototype creation in the context of resource with uri: %s', projection_tree.projection.uri)

        # For every prototype in collection try to create corresponding projections
        for key, prototype in prototypes.items():
            # Set current prototype context to current environment for children node to use
            prototype.context = environment
            # Get context of current node from contexts of parent nodes
            context = prototype.get_context()
            context = context[::-1]

            logger.info('Creating projections for a prototype: %s', prototype)

            # TODO: eval is not safe, consider safer alternative, e.g. JsonPath
            URIs = eval(prototype.uri, locals())

            logger.info('Prototype %s has projections on URIs: %s', prototype, URIs)

            # We get projection URIs based on environment and prototype properties
            # Every URI corresponds to projection object
            for uri in URIs:
                # Get content for a projection
                content = self.driver.get_uri_contents_as_dict(uri)
                name = eval(prototype.name, locals())

                child_projection = ProjectionTree(p_name=name, p_uri=uri)
                # Add newly created projection to projection tree
                projection_tree.add_child(child_projection)
                logger.info('Projection created: %s', child_projection)

                if prototype.type == 'directory':
                    # If prototype has children, continue tree building
                    logger.info('Starting attached prototype projection creation for  prototype: %s with children: %s',
                                prototype, prototype.children)
                    if prototype.children:
                        self.create_projection_tree(prototype.children, child_projection)

                elif prototype.type == 'file':
                    # NOTE: content variable is not accessible during file projection creation!
                    # This is the point where metadata can be extracted from the file
                    #   but file content should be accessed in this case
                    child_projection.projection.type = stat.S_IFREG
                    child_projection.projection.size = 1


