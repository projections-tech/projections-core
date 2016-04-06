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

import argparse
import signal

import psycopg2

__author__ = 'abragin'

import daemon
import logging
from logging import config
import os
import sys
import subprocess
from lockfile import pidlockfile
import shutil

import Pyro4

import yaml

LOCK_FILE = '/var/lock/projections'
PID_LOCK_FILE = '/var/lock/projections.pid'
LOG_FILE = 'projections.log'
# This variable is used to hold reference for daemon object then daemon starts and performed double fork
DAEMON = None

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('projection_daemon')


class ProjectionsDaemon(object):
    logger = logging.getLogger('projection_daemon')

    def __init__(self):
        # This dictionary contains mapping of projectors id`s in database to projections parameters
        self.projections = dict()

        with open('database_connection_config.yaml') as y_f:
            database_connection_parameters = yaml.safe_load(y_f)

        database_host = database_connection_parameters['database_host']
        database_port = database_connection_parameters['database_port']
        user_name = database_connection_parameters['user_name']
        user_password = database_connection_parameters['user_password']

        # Opening connection with database
        self.db_connection = psycopg2.connect(database="projections_database",
                                              user=user_name,
                                              password=user_password,
                                              host=database_host,
                                              port=database_port)

        # Setting connection mode of connection
        self.db_connection.autocommit = False

        # Creating cursor, which will be used to interact with database
        self.cursor = self.db_connection.cursor()

        self.cursor.callproc('projections.daemon_get_projections')

        for row in self.cursor:
            # Dealing with psycopg2 inconsistency here, sometimes cursor returns ("",)
            if len(row) == 1:
                continue

            projection_id, projection_name, mount_point, driver, projector_pid, prototype = row

            # Adding projection
            self.projections[projection_id] = {'projector_subprocess': None,
                                               'mount_point': mount_point,
                                               'prototype_path': prototype,
                                               'driver': driver}

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
            current_projections.append(message_template.format(*row[:-1]))
        return current_projections

    def get_prototypes(self):
        logger.info('List prototypes')
        return 'List prototypes'

    def project(self, projection_name, mount_point, prototype, driver):
        """
        This method starts Projector subprocess which builds projection tree and interacts with FUSE
        :param projection_name: name of projection str
        :param mount_point: projection mount point path str
        :param prototype: path to projection prototype str
        :param driver: name of projection driver str
        :return: status str
        """
        # Checking if projection is in projections table
        self.cursor.callproc('projections.get_projection_id_by_name', [projection_name])
        # Fetching id of projection by it`s name
        projection_id = self.cursor.fetchone()[0]

        self.cursor.callproc('projections.check_mount_point', [mount_point])
        # Checking if mount point is occupied
        mount_point_occupied = self.cursor.fetchone()[0]
        if projection_id is not None:
            # TODO add force flag to override this behaviour
            return 'Error: projection named "{}" already created!'.format(projection_name)
        elif mount_point_occupied:
            return 'Error: mount point "{}" is occupied, please choose another mount point!'.format(mount_point)
        elif not os.path.exists(mount_point):
            return 'Error: provided mount point "{}" does not exist!'.format(mount_point)
        else:
            logger.info('Creating projection with name: %s on %s using prototype: %s and driver: %s.', projection_name,
                        mount_point, prototype, driver)
            self.cursor.callproc('projections.daemon_add_projection', [projection_name, mount_point, driver, prototype])

            projection_id = self.cursor.fetchone()[0]

            # Starting projector subprocess and passing projection parameters to it
            projector = subprocess.Popen([sys.executable,
                                          'db_projector.py',
                                          '-p_n', projection_name,
                                          '-p_id', str(projection_id),
                                          '-m', mount_point,
                                          '-p', prototype,
                                          '-d', driver], stdout=subprocess.DEVNULL)
            # Adding projector to projectors
            self.projections[projection_id] = {'projector_subprocess': projector,
                                               'mount_point': mount_point,
                                               'prototype_path': prototype,
                                               'driver': driver}

            logger.debug('Projector PID: %s', projector.pid)
            # Registering projector in database
            self.cursor.callproc("projections.daemon_set_projection_projector_pid", [projection_id, projector.pid])
            self.db_connection.commit()
            return 'Projection "{}" created and started!'.format(projection_name)

    def start(self, projection_name):
        """
        This method starts previously stopped projection which projector where set to None by creation of new projector
        subprocess
        :param projection_name: name of projection to start str
        :return: status str
        """
        logger.info('Starting projection with name: %s', projection_name)
        self.cursor.callproc('projections.get_projection_id_by_name', [projection_name])
        # Fetching id of projection by it`s name
        projection_id = self.cursor.fetchone()[0]
        # Checking if projection exists
        if projection_id is not None:
            # Setting data required to run projection
            projection_data = self.projections[projection_id]
            # If projection`s mount point does not exist skip projections start
            if not projection_data['mount_point']:
                return 'Error: projection "{}" mount point "{}" does not exist!'.format(projection_name,
                                                                                        projection_data['mount_point'])
            # Checking if projection is already managed
            if projection_data['projector_subprocess'] is None:
                projector = subprocess.Popen([sys.executable,
                                              'db_projector.py',
                                              '-p_n', projection_name,
                                              '-p_id', str(projection_id),
                                              '-m', projection_data['mount_point'],
                                              '-p', projection_data['prototype_path'],
                                              '-d', projection_data['driver'],
                                              '-r'],
                                             stdout=subprocess.DEVNULL)
                # Setting projections projector
                self.projections[projection_id]['projector_subprocess'] = projector

                # Registering projector in database
                self.cursor.callproc("projections.daemon_set_projection_projector_pid", [projection_id, projector.pid])
                self.db_connection.commit()

                return 'Projection {} started!'.format(projection_name)
            else:
                return 'Projection {} is already running!'.format(projection_name)
        else:
            return 'Projection {} does not exist!'.format(projection_name)

    def stop(self, projection_name):
        """
        This method stops projection, terminating it`s projector subprocess
        :param projection_name:
        :return: status str
        """
        logger.info('Stopping projection with name: %s', projection_name)
        self.cursor.callproc('projections.get_projection_id_by_name', [projection_name])
        # Fetching id of projection by it`s name
        projection_id = self.cursor.fetchone()[0]

        if projection_id is not None:
            # Checking if projection  is already stopped
            projection_data = self.projections[projection_id]
            if projection_data['projector_subprocess'] is not None:
                # Terminating projector process and setting it to None for projection with this id
                projector_to_stop = self.projections[projection_id]['projector_subprocess']
                projector_to_stop.terminate()
                projector_to_stop.communicate()

                self.projections[projection_id]['projector_subprocess'] = None

                # Updating projection projector`s pid accordingly
                self.cursor.callproc('projections.daemon_stop_projection', [projection_name])
                self.db_connection.commit()
                return 'Stopped projection: "{}".'.format(projection_name)
            else:
                return 'Projection is already stopped!'
        else:
            return 'There is no projection named "{}"'.format(projection_name)

    def remove_projection(self, projection_name):
        """
        This method perform projection removal from database
        :param projection_name: projection to remove name str
        :return: status str
        """
        logger.info('Removing projection with name: %s', projection_name)
        self.cursor.callproc('projections.get_projection_id_by_name', [projection_name])
        # Fetching id of projection by it`s name
        projection_id = self.cursor.fetchone()[0]

        if projection_id is not None:
            if self.projections[projection_id]['projector_subprocess'] is not None:
                # TODO add flag to force remove projection
                logger.debug('Projection %s is currently running, stopping it!', projection_name)
                self.stop(projection_name)

                # Removing fetched projections data
                projection_data_dir = os.path.join('.projections_contents', projection_name)
                if os.path.exists(projection_data_dir):
                    shutil.rmtree(projection_data_dir)

            self.cursor.callproc('projections.daemon_remove_projection', [projection_name])
            self.db_connection.commit()

            return 'Projection {} removed.'.format(projection_name)
        else:
            return 'Error: Attempting to remove non existant projection!'

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
        This method performs search on projection using SQL as query  language
        :param projection_name: name of projection in which to search str
        :param path: path to search str
        :param query: query str
        :return: status str
        """
        logger.info('Perform search on projection: %s on path: %s', projection_name, path)
        try:
            self.cursor.callproc('projections.search_projections', [projection_name, path, query])
            return [row[0] for row in self.cursor]
        except psycopg2.ProgrammingError as e:
            # Perform rollback of database or else postgresql hangs
            self.db_connection.rollback()
            return 'Error: query is malformed: {}'.format(e)

    def stop_daemon(self):
        """
        This method stops current daemon projections and closes daemon database connections
        :return: None
        """
        logger.debug('Shutting down projections and flushing database connections!')

        self.cursor.callproc('projections.daemon_get_projections')

        for row in self.cursor:
            # Dealing with psycopg2 inconsistency here, sometimes cursor returns ("",)
            if len(row) == 1:
                continue

            logger.debug('Projection to stop at shutdown: %s', row)
            projection_id, projection_name, mount_point, driver, projector_pid, prototype = row
            if projector_pid is None:
                continue
            self.stop(projection_name)

        self.cursor.close()
        self.db_connection.close()


def start_daemon_listener():
    """
    This function creates projections daemon object and registers it in PYRO daemon, which starts requests loop
    :return: None
    """
    global DAEMON
    # This function will be called as child process so we need to initialize loggers again
    logging.config.fileConfig('logging.cfg')
    logger = logging.getLogger('projection_daemon')

    pyro_daemon = Pyro4.Daemon()

    projections_daemon = ProjectionsDaemon()
    # Setting DAEMON global variable, which will be used to perform clean up in stop_daemon function
    DAEMON = projections_daemon

    uri = pyro_daemon.register(projections_daemon)

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
    global DAEMON
    # Perform daemon clean up, stopping running projections and closing daemon database connection
    DAEMON.stop_daemon()
    sys.exit()


# To start projection daemon simply type: ./projections_daemon.py -start
# Then use projections_cli.py client to send command to running daemon.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Projection Daemon')

    options = parser.add_mutually_exclusive_group(required=True)
    options.add_argument('-start', action='store_true', help='Start projections daemon.')
    options.add_argument('-stop', action='store_true', help='Stop projections daemon.')

    args = parser.parse_args()

    if args.start:
        if not os.path.exists(PID_LOCK_FILE):
            # Creating pid lock file which will be used by daemon context manager.
            # This lock file is removed when daemon stops.
            lock_file = pidlockfile.PIDLockFile(PID_LOCK_FILE)
            logger.info('Starting projections daemon!')
            # Create context. For documentation see: https://www.python.org/dev/peps/pep-3143/
            context = daemon.DaemonContext(
                    pidfile=lock_file,
                    stdout=open(LOG_FILE, 'wb'),
                    stderr=sys.stdout,
                    working_directory=os.getcwd())
            # Define signal-action mapping
            context.signal_map = {
                signal.SIGTERM: stop_daemon
            }
            # Starting daemon
            with context:
                start_daemon_listener()
        else:
            # If pid lock file exists daemon is already running
            logger.info('Projections daemon is already running!')
    elif args.stop:
        # If pid lock file exists daemon is running and we can send termination signal to it
        if os.path.exists(PID_LOCK_FILE):
            logger.info('Stoppping projections daemon!')
            with open(PID_LOCK_FILE) as pid_file:
                daemon_pid = pid_file.readline()
            os.kill(int(daemon_pid), 15)
            logger.info('Projections daemon is stopped!')
        else:
            logger.info('Projections daemon is not running!')
