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

from argparse import ArgumentParser

__author__ = 'abragin'

import Pyro4
from Pyro4.errors import CommunicationError

#
# Command-line interface for interaction with projecitons daemon.
#
# Main task of this script is to parse arguments and call projection daemon with appropriate commands.
#

# This file is used for storage of PYRO URI.
# This behavior may be changed in future releases due to multiuser usage patterns.
LOCK_FILE = '/var/lock/projections'

# Maps CLI commands to daemon methods.
# Simplifies API by removing command aliases. If no mapping is provided superuser is used as a method name.
COMMAND_METHOD_MAPPING = {
    'ps': 'get_projections',
    'projections': 'get_projections',
    'prototypes': 'get_prototypes',
    'pt': 'get_prototypes',
    'drivers': 'get_drivers',
    'rm': 'remove_projection',
    'rmp': 'remove_prototype'
}


def _create_args_parser():
    """
    Create main parser and add subparsers to it. Each subparser corresponds to command.

    :return: parser object to use in main method
    """
    parser = ArgumentParser(description='Call Projections meta filesystem commands.')
    subparsers = parser.add_subparsers(title='Commands for Projections operation',
                                       description='valid commands for interaction with Projections objects',
                                       help='To view detailed command description type <command> -h',
                                       dest='command',
                                       metavar='')

    parser_ps = subparsers.add_parser('projections', aliases=['ps'], help='List projections.')

    parser_prototypes = subparsers.add_parser('prototypes', aliases=['pt'], help='List prototypes.')

    parser_project = subparsers.add_parser('project', help='Create new projection.')
    # Random name generation is not recommended since projection names must be unique,
    # so there is a chance of name clashed if generated outside database.
    parser_project.add_argument('-n', '--projection_name', required=True,
                                help='Name of a projection to be created (required).')
    parser_project.add_argument('-m', '--mount_point', required=True,
                                help='Folder in a system to mount projection to. '
                                     'Should exists and be empty (the latter is subject to change in future releases).')
    parser_project.add_argument('-p', '--prototype', required=True,
                                help='Path to prototype file to create projection for.')
    parser_project.add_argument('-d', '--driver', required=True,
                                help='Name of the driver to use with projection')

    parser_start = subparsers.add_parser('start', help='Start projection.')
    parser_start.add_argument('-n', '--projection_name', required=True,
                              help='Name of a projection to start (required).')

    parser_stop = subparsers.add_parser('stop', help='Stop projection.')
    parser_stop.add_argument('-n', '--projection_name', required=True,
                             help='Name of a projection to stop (required).')

    parser_rm = subparsers.add_parser('remove_projection', aliases=['rm'], help='Remove projection.')
    parser_rm.add_argument('-n', '--projection_name', required=True,
                           help='Name of projection to remove.')

    parser_rmp = subparsers.add_parser('remove_prototype', aliases=['rmp'], help='Remove projection prototype.')
    parser_rmp.add_argument('-n', '--projection_name', required=True,
                            help='Name of prototype to remove.')

    parser_drivers = subparsers.add_parser('drivers',
                                           help='List drivers for resources projection registered in the system.')

    parser_search = subparsers.add_parser('search', help='Search projections metadata.')
    parser_search.add_argument('-p', '--path', required=True, help='Limit search to specified path.')
    parser_search.add_argument('-q', '--query', required=True, help='Search objects conforming to the query.')
    parser_search.add_argument('-n', '--projection_name', required=True,
                               help='Name of a projection where to search.')

    return parser


def _get_daemon_object():
    try:
        with open(LOCK_FILE, 'r') as f:
            uri = f.readline()
            projections_daemon = Pyro4.Proxy(uri)
            return projections_daemon
    except Exception as e:
        print('No Projections daemon could be located in the system. Command will be terminated.')
        return None


if __name__ == '__main__':
    # Check whether projection daemon could be accessed
    projection_daemon = _get_daemon_object()

    # Create parser and check arguments
    parser = _create_args_parser()
    args = parser.parse_args()
    if not args.command:
        parser.print_help()

    # Get command and command arguments from argument parser
    command = args.command
    command_arguments = vars(args)

    # Replace commands with method names
    if command in COMMAND_METHOD_MAPPING:
        command = COMMAND_METHOD_MAPPING[command]

    # Remove 'command' from arguments dictionary
    del command_arguments['command']

    # print('Calling command {} with arguments: {}'.format(command, command_arguments))

    # This is DEMO implementation that greatly reduces boilerplate code
    # It requires client calls to be strictly compatible with daemon methods.
    # Extra checks (e.g. syntactical checks) may be incorporated between client commands and daemon methods.
    if command is not None:
        try:
            daemon_command = getattr(projection_daemon, command)
            result = daemon_command(**command_arguments)
            if not isinstance(result, list):
                result = [result]
            for res in result:
                print(res)
        except CommunicationError as e:
            print('FATAL: Error communicating with Projections daemon. Make sure that it is up and running! Error: %s',
                  e)
        except Exception as e:
            print('FATAL: Error occurred while trying to connect to daemon:', e)
