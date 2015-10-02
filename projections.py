__author__ = 'abragin'

import io
import logging
import os
from os import path
import stat
import time
import threading

logger = logging.getLogger('projections')

class Projection(object):
    """
    Class for holding projection data. May be subclassed for area-specific projection implementations
    """

    #TODO: consider composite objects here (such as directories). Need to define schema for this.

    def __init__(self, path, uri):
        """
        Create projection object.

        This object is logically immutable. The default variant is 'directory' projection

        :param path: specify path relative to mount point. Path is always prefixed with root '/' sign
        :param uri: resource identifier projection is pointing to
        :return:
        """
        # TODO: replace 'path' semantics with 'name' semantic
        self.path = path
        _, self.name = os.path.split(self.path)
        self.uri = uri
        # Set to nonzero to trigger projection reading
        self.size = 4096
        # File or folder type. Default is file
        self.type = stat.S_IFDIR

    def __str__(self):
        return 'Projection from {} to {}'.format(self.uri, self.path)


class ProjectionTree(object):
    """
    Class for holding tree-like projection structure.
    """

    def __init__(self):
        # TODO: the self.projections and self.children dictionaries should be synchronized!
        # TODO: consider elements removal. It seems that 'true' tree-like structure is needed here
        # Maps paths to projections in projection tree
        self.projections = {}
        # Element for holding projections children
        #self.children = {'/' : set()}
        # Element for holding root element of projection tree
        self.root = None

        self.lock = threading.Lock()

    def get_projection(self, path):
        if path in self.projections:
            return self.projections[path]
        else:
            return None

    def get_children(self, path):
        if path not in self.children:
            raise KeyError('Path is not managed by projection: {}'.format(path))
        return [p.name for p in self.children[path]]

    def add_projection(self, projection, parent):
        self.lock.acquire()
        self.projections[projection.path] = projection
        #self.children[parent.path].add(projection)
        self.lock.release()

class ProjectionPrototype:

    def __init__(self, type):
        self.children = {}
        self.type = type
        self.name = None
        self.uri = None

    def __str__(self):
        return 'ProjectionPrototype for {}'.format(self.uri)

    def __repr__(self):
        return self.__str__()

class Projector:

    def __init__(self, driver):
        # Resource driver
        self.driver = driver

    def create_projection_tree(self, prototypes, projection_tree, parent_projection=None):
        """
        Creates projection tree for a given collection of prototypes
        """
        logger.info('Creating projection tree with a prototypes: %s starting from: %s', prototypes, parent_projection)

        environment = None
        content = None

        if not parent_projection:
            # Create root projection
            logger.info('Creating root node for projection.')
            projection = Projection('/', None)
            projection.size = 4096
            projection.type = stat.S_IFDIR

            self.create_projection_tree(prototypes, projection_tree, projection)
            logging.info('Projection tree creation done.')
            return projection_tree

        # Get resource content, it should be iterable of JSON objects
        # This is environment in which projections are created
        # In many cases it means double request to parent projection resource so it should be optimized
        environment = self.driver.get_content(parent_projection.uri)

        logger.info('Starting prototype creation in the context of resource with uri: %s', parent_projection.uri)

        for prototype in prototypes:
            logger.info('Creating projections for a prototype: %s', prototype)
            # For every prototype in collection try to create corresponding projections
            URIs = eval(prototype.uri)
            logger.info('Prototype %s has projections on URIs: %s', prototype, URIs)

            for uri in URIs:

                # Get content for a projection
                content = self.driver.get_content(uri)

                logger.info('ENV: %s, CONTENT: %s', environment, content)
                logger.info('Name: %s', prototype.name)
                name = eval(prototype.name)

                projection_path = os.path.join(parent_projection.path, name)
                projection = Projection(projection_path, uri)

                # TODO: Consider safer alternative, e.g. JsonPath
                # Prototype name evaluation is based on content object search
                projection.name = name

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
                    # name = eval(prototype.name)
                    #
                    # projection_path = os.path.join(parent_projection.path, name)
                    #
                    # projection = Projection(projection_path, uri)
                    #
                    # # TODO: Consider safer alternative, e.g. JsonPath
                    # # Prototype name evaluation is based on content object search
                    # projection.name = name
                    #
                    # projection_tree.add_projection(projection, parent_projection)
                    #
                    # logger.info('Projection created: %s', projection)

                    projection.type = stat.S_IFREG
                    projection.size = 1




class ProjectionManager(object):
    """
    Should be subclassed to provide real-life implementations.
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
