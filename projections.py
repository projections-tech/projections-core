__author__ = 'abragin'

import logging
import stat
import time

class Projection(object):
    """
    Class for holding projection data. May be subclassed for area-specific projection implementations
    """

    def __init__(self, path, uri):
        """
        Create projection object.

        :param path: specify path on mount point. Path is always prefixed with root '/' sign
        :param uri: resource identifier projection is pointing to
        :return:
        """
        self.path = path
        self.uri = uri

    def __str__(self):
        return 'Projection from {} to {}'.format(self.uri, self.path)

    def get_content(self):
        return bytes('Hello world!')


class ProjectionManager(object):

    def __init__(self):
        logging.info('Creating projection manager')
        self.projections = dict()
        self.resources = dict()

        projections = self.create_projections()

        for p in projections:
            self.projections[p.path] = p
            self.resources[p.uri] = p

        logging.info('Projections: %s, resources: %s', self.projections, self.resources)

    def create_projections(self):
        """
        This method should be overriden in implementation specific manner

        :return: list of Projection objects
        """
        # This list is filled up by projection specific manner
        projections = []

        # TODO: replace stub with real code here
        projection = Projection('/projection', '/uri:parseq.pro/hello')
        projections.append(projection)

        return projections

    def is_managing_path(self, path):
        return path in self.projections or path in self.resources

    def get_resource(self, uri):
        """
        This method should be overriden in implementation specific manner.

        :param uri: uri to get resource from.
        :return: resource content
        """
        logging.info('Requesting resource for uri: %s', uri)
        return b'Hello world!\n'

    def get_projections(self):
        logging.info('Requesting for projections')
        return self.projections.keys()

    def get_link(self, path):
        if path in self.projections:
            return self.projections[path].uri
        else:
            return None

    def get_attributes(self, path):
        if path in self.projections:
            return self.get_projection_attributes(self.projections[path])
        elif path in self.resources:
            return self.get_resource_attributes(self.resources[path])
        else:
            return None

    def get_projection_attributes(self, projection):

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
        # Set type to link anf grant full access to everyone
        # TODO: Check why links can't be readed
        attributes['st_mode'] = (stat.S_IFREG | 0o0777)
        # Set number of hard links to 0
        attributes['st_nlink'] = 0
        # Set id as inode number.
        attributes['st_ino'] = 1

        return attributes

    def get_resource_attributes(self, projection):

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
