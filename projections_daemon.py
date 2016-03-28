#!/usr/bin/env python3

# Copyright 2016  Anton Bragin, Victor Svekolkin
#
# This file is part of Projections.
#
# Projections is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Projections is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Projections.  If not, see <http://www.gnu.org/licenses/>.
import signal

__author__ = 'abragin'

import daemon
import logging
from logging import config
import os
import sys
import time

import Pyro4

LOCK_FILE = '/var/lock/projections'
LOG_FILE = 'projections.log'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('projection_daemon')

class ProjectionsDaemon(object):

    logger = logging.getLogger('projection_daemon')

    # TODO: substitute with real implementation!

    def get_projections(self):
        logger.info('List projections')
        return 'List projections'

    def get_prototypes(self):
        logger.info('List prototypes')
        return 'List prototypes'

    def project(self, projection_name, mount_point, prototype, driver):
        logger.info('Creating projection with name: %s on %s using prototype: %s and driver: %s.')
        return 'Creating projection'

    def start(self, projection_name):
        logger.info('Starting projection with name: %s')
        return 'Starting projection'

    def stop(self, projection_name):
        logger.info('Stopping projection with name: %s')
        return 'Stopping projection'

    def remove_projection(self, projection_name):
        logger.info('Removing projection with name: %s')
        return 'Removing projection'

    def remove_prototype(self, prototype_name):
        logger.info('Removing prototype with name: %s')
        return 'Removing prototype'

    def get_drivers(self):
        logger.info('List drivers')
        return 'List drivers'

    def search(self, path, query):
        logger.info('Do search')
        return 'searching on path: {} with query: {}'.format(path, query)


def start_daemon_listener():
    # This function will be called as child process so we need to initialize loggers again
    logging.config.fileConfig('logging.cfg')
    logger = logging.getLogger('projection_daemon')

    pyro_daemon = Pyro4.Daemon()
    uri = pyro_daemon.register(ProjectionsDaemon)

    with open(LOCK_FILE, 'w') as f:
        f.write(str(uri))
        f.write('\n')

    logger.info('Projections daemon is starting. URI: {}'.format(uri))
    pyro_daemon.requestLoop()

def stop_daemon(signum, frame):
    """
    This method is called when Projections daemon receives request to terminate.
    The signature of method is fixed. See: https://docs.python.org/2/library/signal.html#signal.signal

    :param signum: signal number
    :param frame: current stack frame
    """
    logger.info('Signal to stop daemon received. Terminating projections daemon.')

    # TODO: implement all cleaning up stuff here.
    # Terminate projections FUSE subprocesses, flush database connections, etc.


# To start projection daemon simple type: ./projections_daemon.py
# Then use projections_cli.py client to send command to running daemon.
if __name__ == '__main__':
    logger.debug('Starting daemon')
    # Create context. For documentation see: https://www.python.org/dev/peps/pep-3143/
    context = daemon.DaemonContext(
            pidfile=open('/var/lock/projections.pid', 'wb'),
            stdout=open(LOG_FILE, 'wb'),
            stderr=sys.stdout,
            working_directory=os.getcwd())

    context.signal_map = {
        signal.SIGTERM: stop_daemon
    }

    with context:
        start_daemon_listener()
