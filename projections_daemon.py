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

import getpass
import signal

import psycopg2

__author__ = 'abragin'

import daemon
import logging
from logging import config
import os
import sys
import subprocess

import Pyro4

LOCK_FILE = '/var/lock/projections'
LOG_FILE = 'projections.log'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('projection_daemon')


class ProjectionsDaemon(object):
    logger = logging.getLogger('projection_daemon')

    def __init__(self):
        # This dictionary contains mapping of projectors id`s in database to running subprocessess objects
        self.projectors = dict()
        self.logger.debug('Id of projectors dict at init: %s', id(self.projectors))
        # Opening connection with database
        self.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))

        # Setting connection mode of connection
        self.db_connection.autocommit = False

        # Creating cursor, which will be used to interact with database
        self.cursor = self.db_connection.cursor()

    def get_projections(self):
        """
        This method returns a list of projections in database
        :return: list of projection names strings
        """
        self.cursor.callproc('projections.daemon_get_projections')

        current_projections = list()
        message_template = 'Projection id: {0} Projection name: {1} Mount point: {2} Driver: {3} Projector PID: {4}'

        # Format message template and return projections list
        for row in self.cursor:
            current_projections.append(message_template.format(*row))
        return current_projections

    def get_prototypes(self):
        logger.info('List prototypes')
        return 'List prototypes'

    def project(self, projection_name, mount_point, prototype, driver):
        """

        :param projection_name:
        :param mount_point:
        :param prototype:
        :param driver:
        :return:
        """
        # Checking if projection is in projections table
        self.cursor.callproc('projections.get_projection_id_by_name', [projection_name])
        # Fetching id of projection by it`s name
        fetch_result = self.cursor.fetchone()

        if fetch_result is not None:
            # TODO add force flag to override this behaviour
            return 'Projection already created!'
        else:
            logger.info('Creating projection with name: %s on %s using prototype: %s and driver: %s.', projection_name,
                        mount_point, prototype, driver)
            self.cursor.callproc('projections.daemon_add_projection', [projection_name, mount_point, driver])

            projection_id = self.cursor.fetchone()[0]

            projector = subprocess.Popen([sys.executable,
                                          'db_projector.py',
                                          '-p_id', str(projection_id),
                                          '-m', mount_point,
                                          '-p', prototype,
                                          '-d', driver], stdout=subprocess.DEVNULL)

            self.projectors[projection_id] = {'projector_subprocess': projector,
                                              'mount_point': mount_point,
                                              'prototype_path': prototype,
                                              'driver': driver}

            self.db_connection.commit()

        logger.debug('Projector PID: %s', projector.pid)

        self.cursor.callproc("projections.daemon_set_projection_projector_pid", [projection_id, projector.pid])
        self.logger.debug('Id of projectors dict at project: %s', id(self.projectors))
        self.db_connection.commit()
        return 'Projection {} created and started!'.format(projection_name)

    def start(self, projection_name):
        """

        :param projection_name:
        :return:
        """
        logger.info('Starting projection with name: %s', projection_name)
        self.cursor.callproc('projections.get_projection_id_by_name', [projection_name])
        # Fetching id of projection by it`s name
        fetch_result = self.cursor.fetchone()
        # Checking if projection exists
        if fetch_result is not None:
            projection_id = fetch_result[0]
            # Setting data required to run projection
            projection_data = self.projectors[projection_id]

            if projection_data['projector_subprocess'] is None:
                projector = subprocess.Popen([sys.executable,
                                              'db_projector.py',
                                              '-p_id', str(projection_id),
                                              '-m', projection_data['mount_point'],
                                              '-p', projection_data['prototype_path'],
                                              '-d', projection_data['driver'],
                                              '-r'],
                                             stdout=subprocess.DEVNULL)

                self.projectors[projection_id]['projector_subprocess'] = projector
                return 'Projection {} started!'.format(projection_name)
            else:
                return 'Projection {} is already running!'.format(projection_name)
        else:
            return 'Projection {} does not exist!'.format(projection_name)

    def stop(self, projection_name):
        """

        :param projection_name:
        :return:
        """
        logger.info('Stopping projection with name: %s', projection_name)
        self.cursor.callproc('projections.get_projection_id_by_name', [projection_name])
        # Fetching id of projection by it`s name
        fetch_result = self.cursor.fetchone()

        if fetch_result is not None:
            projection_id = fetch_result[0]
            # Terminating projector process and setting it to None for projection with this id
            projector_to_stop = self.projectors[projection_id]['projector_subprocess']
            projector_to_stop.terminate()
            self.projectors[projection_id]['projector_subprocess'] = None

            # Updating projection projector`s pid accordingly
            self.cursor.callproc('projections.daemon_stop_projection', [projection_name])
            self.db_connection.commit()
            return 'Stopped projection: {}'.format(projection_name)
        else:
            return 'There is no projection named {}'.format(projection_name)

    def remove_projection(self, projection_name):
        """

        :param projection_name:
        :return:
        """
        logger.info('Removing projection with name: %s', projection_name)
        self.cursor.callproc('projections.get_projection_id_by_name', [projection_name])
        # Fetching id of projection by it`s name
        projection_id = self.cursor.fetchone()[0]

        if self.projectors[projection_id]['projector_subprocess'] is not None:
            # TODO add flag to force remove projection
            logger.debug('Projection %s is currently running, stopping it!', projection_name)
            self.stop(projection_name)

        self.cursor.callproc('projections.daemon_remove_projection', [projection_name])
        return 'Projection {} removed.'.format(projection_name)

    def remove_prototype(self, prototype_name):
        """

        :param prototype_name:
        :return:
        """
        logger.info('Removing prototype with name: %s')
        return 'Removing prototype'

    def get_drivers(self):
        """

        :return:
        """
        logger.info('List drivers')
        return 'List drivers'

    def search(self, projection_name, path, query):
        """

        :param projection_name:
        :param path:
        :param query:
        :return:
        """
        logger.info('Do search')
        self.cursor.callproc('projections.search_projections', [projection_name, path, query])
        return [row[0] for row in self.cursor]

    def stop_daemon(self):
        """

        :return:
        """
        for projection_id, projection_data in self.projectors.items():
            projector = projection_data['projector_subprocess']

            if projector is None:
                continue
            logger.debug('Terminating projection: %s, projector: %s', projection_id, projector)
            projector.terminate()

        self.cursor.close()
        self.db_connection.close()


def start_daemon_listener():
    """

    :return:
    """

    # This function will be called as child process so we need to initialize loggers again
    logging.config.fileConfig('logging.cfg')
    logger = logging.getLogger('projection_daemon')

    pyro_daemon = Pyro4.Daemon()

    uri = pyro_daemon.register(ProjectionsDaemon())

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
