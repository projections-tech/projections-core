__author__ = 'abragin'

import io
import logging
import os
import stat
import time
import threading
import yaml

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('projections')


class Tree(object):
    """
    Generic tree class which holds common methods for "tree" structure, must be subclassed to create Prototype and
    Projection classes
    """
    def __init__(self, name=None, data=None):
        """
        Initializes Tree node with optional arguments name and data
        :param name: Name of current node
        :param data: Field that holds data associated with current node
        """
        self.name = name
        self.data = data
        self.parent = None
        self.children = {}

    def __split_path_in_parts(self, path):
        """
        Split path string in parts preserving path root
        :param path: path string
        :return: list of path parts strings
        """
        temp = []
        while True:
            head, tail = os.path.split(path)
            if head == '/' or head == '':
                temp.append(tail)
                temp.append(head)
                return temp[::-1]
            temp.append(tail)
            path = head

    def add_child(self, node):
        node.parent = self
        self.children[node.name] = node

    def tree_to_list(self):
        """
        Converts list of lists tree to Tree of Trees
        """
        return [self.data]+[c.tree_to_list() for c in self.children.values()]

    def find(self, name):
        """
        Find node in tree according to its name
        :param value:
        :return: node in Tree, which is Tree object, None if no corresponding node where found
        """

        if self.name == name:
            return self
        else:
            for c in self.get_children():
                result = c.find(name)
                if not result is None:
                    return result
            return None

    def get_path(self):
        """
        Returns nodes list from root to current node.
        """
        result = []
        parent, child = self.parent, self
        while parent:
            result.append(parent)
            parent, child = parent.parent, parent
        return result[::-1]

    def get_children(self):
        """
        Returns sorted by name list of children of current node
        """
        return sorted(list(self.children.values()), key=lambda x: x.name)

    def traverse_pre_order(self):
        """
        Used to traverse tree in pre order manner
        """
        yield self
        for c in self.get_children():
            for v in c.traverse_pre_order():
                yield v

    def traverse_post_order(self):
        """
        Used to traverse tree in post order manner
        """
        for c in self.get_children():
            for v in c.traverse_post_order():
                yield v
        yield self

    def traverse_depth_first(self):
        """
        Used to traverse tree depth first
        """
        to_yield = [self]
        while to_yield:
            node = to_yield.pop()
            for c in node.get_children():
                to_yield.append(c)
            yield node

    def find_node_by_path(self, path_to_node):
        """
        Finds node in a tree according to path, by iterating path parts on part at a time on level of tree structure,
        per iteration, assuming that starting node is root of path_to_node.
        """
        path = self.__split_path_in_parts(path_to_node)
        temp_node = self
        if temp_node.name == path[0]:
            for item in path[1:]:
                if item in temp_node.children:
                    temp_node = temp_node.children[item]
                else:
                    return None
        else:
            return None
        return temp_node


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

    This is class STUB.

    TODO: Implementation that allows path traversal and children retrieval should be created.
    Class functionality is highly similar to ProjectionManager behavior, the latter can achieve this by composition.

    """

    def __init__(self):
        # TODO: consider elements removal. It seems that 'true' tree-like structure is needed here
        # Maps paths to projections in projection tree
        self.projections = {}
        self.root = None

        # Be aware of atomicity of operations
        self.lock = threading.Lock()

    def get_projection(self, path):
        """
        Get projection for a path provided.

        :param path: path the projection resides on
        :return: projection or None
        """
        if path in self.projections:
            return self.projections[path]
        else:
            return None

    def get_children(self, path):
        """
        Return first-order descendants of the projection on the path specified.
        Implements 'get-directory-content-like' functionality.

        :param path: path to retrieve child projections for
        :return: list of child projections that can be empty
        """
        raise NotImplemented('This method is not currently implemented')

    def add_projection(self, projection, parent):
        """
        Adds projection to the current tree.

        :param projection: projection object to add.
        :param parent: parent projection the object should be attached to.
        """
        # This should be done atomically regardless of implementation details.
        self.lock.acquire()
        self.projections[projection.path] = projection

        self.lock.release()


class ProjectionPrototype(Tree):
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
        # Initialize Tree class, passing current object as Tree data field
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
        return [node.context for node in self.get_path()]


class PrototypeDeserializer(object):
    """
    Used to deserialize Prototype tree from YAML configuration files
    """
    def __init__(self, data_path):
        """
        Initialize class, passing path to YAML configuration file
        """

        self.prototype_tree, self.resource_uri = self.read_projections(data_path)

    def get_prototypes_tree(self, yaml_dict, parent=None):
        """
        Constructs ProjectionPrototype tree from root dictionary
        :param yaml_dict: dict for root element
        :param parent: parent prototype for tree node
        :return root ProjectionPrototype object
        """
        t = ProjectionPrototype(type=yaml_dict['type'])
        t.name = yaml_dict['name']
        t.uri = yaml_dict['uri']
        t.parent = parent
        if isinstance(yaml_dict['children'], dict):
            t.children = {x[0]: self.get_prototypes_tree(x[1], parent=t) for x in yaml_dict['children'].items()}
            return t
        else:
            return t

    def read_projections(self, data_path):
        """
        Read YAML configuration file on data_path and construct Prototype tree from it
        :param data_path
        :return root ProjectionPrototype object
        """
        with open(data_path) as y_f:
            yaml_dict = yaml.safe_load(y_f)
        return self.get_prototypes_tree(yaml_dict['root']), yaml_dict['resource_uri']


class ProjectionDriver(object):
    """
    Object that has get_content(uri) method returning array of Python dictionaries.
    """

    def get_content(self, uri):
        raise NotImplemented('Implement data retrieval from some projection backend.')


class Projector:
    """
    Class that creates Projection object and assembles them in ProjectionTree object using Prototypes object.
    """

    def __init__(self, driver):
        """
        Create projector object.

        :param driver: driver object that has get_content(uri) method returning array of Python dictionaries
            that are used for ProjectionTree composition.
        """
        assert isinstance(driver, ProjectionDriver), 'Check that driver object is subclass of ProjectionDriver'
        self.driver = driver

    def fetch_context(self, uri):
        """
        Used to fetch contents of uri in microcode
        :param uri: URI
        :return: uri contents
        """
        return self.driver.get_content(uri)

    def create_projection_tree(self, prototypes, projection_tree, parent_projection=None):
        """
        Creates projection tree for a given collection of prototypes.

        Method signature is influenced by recursion implementation and is subject to change.

        :param prototypes: collection of prototypes that are attached to parent projection's prototype
        :param projection_tree: projection tree to extend. Should not be None
        :param parent_projection: parent projection that is root to tree under creation.
            Should exists in projection tree provided.
        :return: void. The provided projection tree object expected to be used.
        """
        logger.info('Creating projection tree with a prototypes: %s starting from: %s', prototypes, parent_projection)

        # Dictionaries that act as a evaluation context for projections. environment is available for both directory and
        #   file prototypes, while content for directory prototypes only
        environment = None
        content = None
        fetch_context = self.fetch_context
        path = os.path

        # This is environment in which projections are created (parent_projection content)
        # TODO: in many cases it means double request to parent projection resource so it should be optimized
        environment = self.driver.get_content(parent_projection.uri)
        logger.debug(environment)
        logger.info('Starting prototype creation in the context of resource with uri: %s', parent_projection.uri)

        # For every prototype in collection try to create corresponding projections
        for key, prototype in prototypes.items():
            # Set current prototype context to current environment for children node to use
            prototype.context = environment
            # Get context of current node from contexts of parent nodes
            context = prototype.get_context()
            context = context[::-1]
            logger.info('Prototype context: %s', len(context))
            logger.info('Creating projections for a prototype: %s', prototype)

            # TODO: eval is not safe, consider safer alternative, e.g. JsonPath
            URIs = eval(prototype.uri, locals())

            logger.info('Prototype %s has projections on URIs: %s', prototype, URIs)

            # We get projection URIs based on environment and prototype properties
            # Every URI corresponds to projection object
            for uri in URIs:
                # Get content for a projection
                content = self.driver.get_content(uri)
                logger.debug('ENV: %s, CONTENT: %s', environment, content)
                name = eval(prototype.name, locals())

                # This may be reconsidered with ProjectionTree implementation
                projection_path = os.path.join(parent_projection.path, name)

                projection = Projection(projection_path, uri)
                projection.name = name

                # Add newly created projection to projection tree
                projection_tree.add_projection(projection, parent_projection)
                logger.info('Projection created: %s', projection)

                if prototype.type == 'directory':

                    projection.type = stat.S_IFDIR
                    projection.size = 4096

                    # If prototype has children, continue tree building
                    logger.info('Starting attached prototype projection creation for  prototype: %s with children: %s',
                                prototype, prototype.children)
                    if prototype.children:
                        self.create_projection_tree(prototype.children, projection_tree, projection)

                elif prototype.type == 'file':
                    # NOTE: content variable is not accessible during file projection creation!
                    # This is the point where metadata can be extracted from the file
                    #   but file content should be accessed in this case
                    projection.type = stat.S_IFREG
                    projection.size = 1


class ProjectionManager(object):
    """
    Should be subclassed to provide real-life implementations.

    TODO: ProjectionManager should include Projector functionality and present it to subclasses where applicable.
    """

    def __init__(self):
        logger.info('Creating projection manager')
        self.projections = dict()
        self.resources = dict()

        projections = self.create_projections()

        for p in projections:
            self.projections[p.path] = p
            self.resources[p.uri] = p

        logger.debug('Projections: %s, resources: %s', self.projections, self.resources)

    def create_projections(self):
        """
        This method should be overriden in implementation specific manner

        :return: list of Projection objects
        """
        # This list is filled up by projection specific manner
        projections = []

        # TODO: replace stub with real code here
        projection = Projection('/projection', 'uri:parseq.pro')
        projection.type = stat.S_IFREG
        projection.size = 1
        projections.append(projection)

        return projections

    def is_managing_path(self, path):
        return path in self.projections or path in self.resources

    def get_resource(self, uri):
        """
        This method should be overriden in implementation specific manner.

        :param uri: resource identifier to get resource from.
        :return: resource content
        """
        # TODO: implement resource downloading

        content = b'Hello World!\n'

        logger.info('Requesting resource for uri: %s', uri)
        return content

    def open_resource(self, path):
        # TODO: implement file header operations
        # TODO: implement resource downloading if not already opened (use header to check)!

        content = b'Hello World!\n'
        resource_io = io.BytesIO(content)

        if path in self.projections:
            self.projections[path].size = len(content)

        logger.info('Resource content size for uri: %s set to: %s', path, self.projections[path].size)

        file_handler = 3

        return file_handler, resource_io

    def remove_resource(self, path):
        if path in self.projections:
            self.projections[path].size = 1

    def get_projections(self, path):
        logger.info('Requesting for projections')
        # TODO: this is hardcoded behavior with only one level of projections. Review it!
        if path == '/':
            return self.projections.keys()
        else:
            return []

    def get_link(self, path):
        if path in self.projections:
            return self.projections[path].uri
        else:
            return None

    def get_attributes(self, path):
        logger.info('Request attributes on path: %s. Projections: %s, resources: %s', path, self.projections, self.resources)
        if path in self.projections:
            return self.get_projection_attributes(self.projections[path])
        elif path[1:] in self.resources:
            return self.get_resource_attributes(self.resources[path[1:]])
        else:
            logging.info('No object for attributes found')
            return None

    def get_xattr(self, path, name):
        logger.info('Requestin projections \'%s\' from path: %s', name, path)
        return self.projections[path].xattrs[name]

    def get_projection_attributes(self, projection):
        logger.info('Get projection attributes')
        now = time.time()
        attributes = dict()

        # Set projection attributes

        # This is implementation specific and should be binded to projector data
        attributes['st_atime'] = now
        # This may be implemented as last projection cashing time is casing is enabled
        attributes['st_mtime'] = now
        # On Unix this is time for metedata modification we can use the same conception
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

    def get_resource_attributes(self, projection):
        logger.info('Get resource attributes')
        now = time.time()
        attributes = dict()

        # Set projection attributes

        # This is implementation specific and should be binded to projector data
        attributes['st_atime'] = now
        # This may be implemented as last projection cashing time is casing is enabled
        attributes['st_mtime'] = now
        # On Unix this is time for metedata modification we can use the same conception
        attributes['st_ctime'] = now
        # If this is projection the size is zero
        attributes['st_size'] = 10
        # TODO: this is implementation specific option! Resource may be directory as well
        # Set type to regular file anf grant full access to everyone
        attributes['st_mode'] = (stat.S_IFREG | 0o0777)
        # Set number of hard links to 0
        attributes['st_nlink'] = 0
        # Set id as inode number.
        attributes['st_ino'] = 1

        return attributes
