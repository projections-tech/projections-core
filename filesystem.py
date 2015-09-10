#!/usr/bin/env python3

__author__ = 'abragin'

import logging
import logging.config
import os
import stat
import sys
import time
from fuse import FUSE, FuseOSError, Operations

from projections import Projection, ProjectionManager

logger = logging.getLogger('projection_filesystem')


class ProjectionFilesystem(Operations):
    def __init__(self, mount_point, data_folder):
        """
        Create projection filesystem for given mountpoint

        :param mount_point: point to mount projection filesystem to
        :param data_folder: pass-through folder for projectins dta storage
        :return: None
        """
        self.mount_point = mount_point
        # This variable defines directory with the actual stored data
        # in future 'mount to itself' behavior can be emulated bu mounting original folder to some place known to the program
        # and mounting projection to the original folder
        self.data_path = os.path.abspath(data_folder)
        logger.info('Activating projection on mount point: %s with data path: %s',
                    self.mount_point, self.data_path)

        logger.debug('Content of data directory: %s', os.listdir(self.data_path))
        # Creating projection manager
        self.projection_manager = ProjectionManager()

    # FUSE methods
    def _extend_data_path(self, path):
        if path.startswith('/'):
            path = path[1:]
        return os.path.join(self.data_path, path)

    def getattr(self, path, fh=None):
        absolute_path = self._extend_data_path(path)
        logger.info('Requesting attributes for path: %s, absolute: %s', path, absolute_path)

        # If attributes belong to projection manager than return
        if self.projection_manager.is_managing_path(path):
            attributes = self.projection_manager.get_attributes(path)
            logger.debug('Projection attributes received: %s', attributes)
            if attributes:
                return attributes

        # Otherwise request attributes from data folder
        st = os.lstat(absolute_path)

        if st:
            logger.debug('Data stats received: %s', st)
            return dict((key, getattr(st, key)) for key in ('st_atime', 'st_ctime',
                                                            'st_gid', 'st_mode', 'st_mtime', 'st_nlink', 'st_size',
                                                            'st_uid'))
        # If nobody manages this path than raise exception
        raise RuntimeWarning('Resource on {} have no associated attributes'.format(path))

    def readdir(self, path, fh):
        logger.info('Reading directory entries. Path: %s, fh: %s', path, fh)
        absolute_path = self._extend_data_path(path)

        dirents = []
        # Append data directory entries
        if os.path.isdir(absolute_path):
            dirents.extend(os.listdir(absolute_path))

        logger.debug('Directory entries: %s', dirents)

        def remove_heading_slash(path):
            if path[0] == '/':
                path = path[1:]
            return path

        # Get projections and exclude the ones that are already on the disk
        projections = [remove_heading_slash(p) for p in self.projection_manager.get_projections(path) if
                       remove_heading_slash(p) not in dirents]
        logger.debug('Filtered projections: %s', projections)
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
        logger.info('Reading link on path: %s', path)
        link_path = self.projection_manager.get_link(path)

        if link_path:
            logger.info('Returning path: %s', link_path)
            return link_path
        else:
            # TODO: implement work with data folder links
            logger.info('Managing path %s by data folder', path)
            return ''

    def read(self, path, length, offset, fh):
        # Request content from real data
        data_path = self._extend_data_path(path)
        if os.path.exists(data_path):
            logger.debug('Requesting content from real path. Path: %s, length: %s, offset: %s, fh :%s', data_path, length, offset, fh)
            os.lseek(fh, offset, os.SEEK_SET)
            return os.read(fh, length)

        # Looking fo data from projection manager
        logger.info('Reading file on path: %s. Length: %s, offset: %s, header: %s', path, length, offset, fh)
        if self.projection_manager.is_managing_path(path):
            logger.debug('Reading projection resource on path: %s', path)
            content = self.projection_manager.get_resource(path)
            logger.debug('Got content for uri: %s: %s', path, content)
            return content

    def open(self, path, flags):
        logger.info('Opening file on path: %s with flags: %s', path, flags)
        # If the path is managed by projection manager than it should place original resource on drive before opening
        if self.projection_manager.is_managing_path(path):
            file_header, resource_io = self.projection_manager.open_resource(path)
            logger.debug('Opening resource at path: %s returned header: %s', path, file_header)

            logger.info('Saving resource content to local drive')
            # Create folder if not exists
            dirname = os.path.dirname(self._extend_data_path(path))

            logger.debug('Creating directory if not exists: %s', dirname)
            if not os.path.exists(dirname):
                logger.debug('Creating directory: %s', dirname)
                os.makedirs(dirname)

            data_path = self._extend_data_path(path)

            with open(data_path, 'wb') as f:
                f.write(resource_io.read())

        # Opening real file that was created
        return os.open(self._extend_data_path(path), flags)

    def rmdir(self, path):
        logging.info('Removing node on path: %s', path)

        if self.projection_manager.is_managing_path(path):
            logging.debug('Path is projection, do nothing')
            return None

        full_path = self._extend_data_path(path)
        return os.rmdir(full_path)

    def unlink(self, path):
        logging.info('Unlinking path: %s', path)
        data_path = self._extend_data_path(path)
        if os.path.exists(data_path):
            if self.projection_manager.is_managing_path(path):
                self.projection_manager.remove_resource(path)

            return os.unlink(data_path)

    # def getxattr(self, path, name, position=0):
    #     logging.info('Get xattr \'%s\' from path: %s', name, path)
    #
    #     if self.projection_manager.is_managing_path(path):
    #         logging.debug('Path is managed by projection manager')
    #         return self.projection_manager.get_xattr()


def main(mountpoint, data_folder, foregrount=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    fuse = FUSE(ProjectionFilesystem(mountpoint, data_folder), mountpoint, foreground=foregrount)
    return fuse


if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    # TODO: replace with argparse
    # Get mount point from args
    if len(sys.argv) != 3:
        print('usage: %s <mountpoint> <data folder>' % sys.argv[0])
        exit(1)

    main(sys.argv[1], sys.argv[2])
