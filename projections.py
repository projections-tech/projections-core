__author__ = 'abragin'

import io
import logging
import stat
import time

logger = logging.getLogger('projections')

class Projection(object):
    """
    Class for holding projection data. May be subclassed for area-specific projection implementations
    """

    #TODO: consider composite objects here (such as directories). Need to define schema for this.

    def __init__(self, path, uri):
        """
        Create projection object.

        :param path: specify path on mount point. Path is always prefixed with root '/' sign
        :param uri: resource identifier projection is pointing to
        :return:
        """
        self.path = path
        self.uri = uri
        # TODO: this is demo implementation. There should be more elegant way to do the same, e.g. via special file types
        # Set to nonzero to trigger projection reading
        self.size = 1
        # File or folder type. Default is file
        self.type = stat.S_IFREG
        self.xattrs = {}

        # Other attributes as well as extended attributes may be stored here
        self.download_process = None

    def __str__(self):
        return 'Projection from {} to {} with xattrs: {}'.format(self.uri, self.path, self.xattrs)


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
