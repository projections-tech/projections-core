#!/usr/bin/env python3

__author__ = 'abragin'

import logging
import os
import stat
import sys
import time
from fuse import FUSE, FuseOSError, Operations

from projections import Projection, ProjectionManager

# This variable defines directory with the actual stored data
# in future 'mount to itself' behavior can be emulated bu mounting original folder to some place known to the program
# and mounting projection to the original folder
DATA_DIRECTORY = 'data'


class ProjectionFilesystem(Operations):

    def __init__(self, mount_point):
        self.mount_point = mount_point
        self.data_path = os.path.abspath(DATA_DIRECTORY)
        logging.info('Activating projection on mount point: %s with data path: %s',
                     self.mount_point, self.data_path)

        logging.debug('Content of data directory: %s', os.listdir(self.data_path))
        # Creating projection manager
        self.projection_manager = ProjectionManager()

    # FUSE methods
    def _extend_data_path(self, path):
        if path.startswith('/'):
            path = path[1:]
        return os.path.join(self.data_path, path)

    def getattr(self, path, fh=None):
        absolute_path = self._extend_data_path(path)
        logging.info('Requesting attributes for path: %s, absolute: %s', path, absolute_path)

        # If attributes belong to projection manager than return
        attributes = self.projection_manager.get_attributes(path)
        logging.debug('Projection attributes received: %s', attributes)
        if attributes:
            return attributes

        # Otherwise request attributes from data folder
        st = os.lstat(absolute_path)

        if st:
            logging.debug('Data stats received: %s', st)
            return dict((key, getattr(st, key)) for key in ('st_atime', 'st_ctime',
                     'st_gid', 'st_mode', 'st_mtime', 'st_nlink', 'st_size', 'st_uid'))

        raise RuntimeWarning('Resource on {} have no associated attributes'.format(path))

    def readdir(self, path, fh):
        logging.info('Reading directory entries. Path: %s, fh: %s', path, fh)
        absolute_path = self._extend_data_path(path)

        dirents = []
        # Append data directory entries
        if os.path.isdir(absolute_path):
            dirents.extend(os.listdir(absolute_path))

        logging.debug('Directory entries: %s', dirents)

        def remove_heading_slash(path):
            if path[0] == '/':
                path = path[1:]
            return path

        # Get projections and exclude the ones that are already on the disk
        # TODO: this is hardcoded behavior with only one level of projections. Review it!
        if path == '/':
            projections = [remove_heading_slash(p) for p in self.projection_manager.get_projections() if p not in dirents]

            logging.debug('Filtered projections: %s', projections)
            # Mix-in projections to directory entries
            dirents.extend(projections)

        for d in dirents:
            yield d

    def readlink(self, path):
        """
        Return a string representing the path to which the symbolic link points.

        The result may be either an absolute or relative pathname.

        :param path:
        :return:
        """
        logging.info('Reading link on path: %s', path)
        link_path = self.projection_manager.get_link(path)

        if link_path:
            logging.info('Returning path: %s', link_path)
            return link_path
        else:
            # TODO: implement work with data folder links
            logging.info('Managing path %s by data folder', path)
            return ''

    def read(self, path, length, offset, fh):
        logging.info('Reading file on path: %s. Length: %s, offset: %s, header: %s', path, length, offset, fh)
        if self.projection_manager.is_managing_path(path):
            logging.debug('Reading projection resource on path: %s', path)
            content = self.projection_manager.get_resource(path)
            logging.debug('Got content for uri: %s: %s', path, content)
            return content
        else:
            with open(self._extend_data_path(path), 'rb') as f:
                content = f.read(length)
                logging.debug('Got content for uri: %s: %s', path, content)
                return content

    def open(self, path, flags):
        logging.info('Opening file on path: %s with flags: %s', path, flags)
        if self.projection_manager.is_managing_path(path):
            file_header = self.projection_manager.open_resource(path)
            logging.debug('Opening resource at path: %s returned header: %s', path, file_header)
            return file_header
        else:
            return os.open(self._extend_data_path(path), flags)

    # Projection methods. Will be moved to separate class
    # def get_projections(self):
    #
    #     projections = dict()
    #     projections['/link'] = Projection(1, 'link')
    #
    #     return projections



    # Filesystem methods
    # def access(self, path, mode):
    #     """
    #     Test for access to path.
    #
    #     See the Unix man page access(2) for more information.
    #
    #     :param path:
    #     :param mode:
    #     :return: Return True if access is allowed, False if not.
    #     """
    #     logging.info('access method called with parameters path : %s, mode : %s', path, mode)
    #     return True

    # def chmod(self, path, mode):
    #     logging.info('chmod method called with parameters path : %s, mode : %s', path, str(mode))
    #
    # def chown(self, path, uid, gid):
    #     logging.info('chown method called with parameters path : %s, uid : %s, gid: %s', path, uid, gid)


    # def readlink(self, path):
    #     """
    #     Return a string representing the path to which the symbolic link points.
    #
    #     The result may be either an absolute or relative pathname.
    #
    #     :param path:
    #     :return:
    #     """
    #     return path

    # def mknod(self, path, mode, dev):
    #     """
    #     Create a filesystem node (file, device special file or named pipe) named filename.
    #
    #     :param path:
    #     :param mode: permissions to use and the type of node to be created
    #     :param dev: device defines the newly created device special file
    #     :return:
    #     """
    #     return os.mknod(os.path.join(self.mount_point, path), mode, dev)
    #
    # def rmdir(self, path):
    #     return super.rmdir(path)
    # 
    # def mkdir(self, path, mode):
    #     return super.mkdir(path, mode)



if __name__ == '__main__':

    # Get mount point from args
    if len(sys.argv) != 2:
        print('usage: %s <mountpoint>' % sys.argv[0])
        exit(1)

    logging.basicConfig(format='#%(levelname)-8s %(name)-12s [%(asctime)s]  %(message)s', level=logging.DEBUG)

    # Specify FUSE mount options as **kwargs here. For value options use value=True form
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    fuse = FUSE(ProjectionFilesystem(sys.argv[1]), sys.argv[1], foreground=True, nonempty=True)
