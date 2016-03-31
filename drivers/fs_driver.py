#!/usr/bin/env python3

#    Copyright 2016  Anton Bragin, Victor Svekolkin
#
#    This file is part of Projections.
#
#    Projections is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Projections is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Projections.  If not, see <http://www.gnu.org/licenses/>.

import logging
import logging.config
import os

from projections import ProjectionDriver

logger = logging.getLogger('filesystem_driver')


class FSDriver(ProjectionDriver):
    def __init__(self, uri, driver_config_path):
        self.uri = uri
        self.driver_configuration = self.read_config(driver_config_path)

    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """
        uri = os.path.join(uri)

        # Directory projection returns list of it`s children as metadata
        if os.path.isdir(uri):
            return {'name': os.path.basename(uri),
                    'resource_uri': os.path.abspath(uri),
                    'size': os.path.getsize(uri),
                    'type': 'dir',
                    'children': [
                        self.get_uri_contents_as_dict(os.path.join(os.path.abspath(uri), p)) for p in os.listdir(uri)
                        ],
                    'extension': None}
        else:
            return {'name': os.path.basename(uri),
                    'resource_uri': os.path.abspath(uri),
                    'size': os.path.getsize(uri),
                    'type': 'file',
                    'extension': os.path.splitext(uri)[1]}

    def get_uri_contents_as_bytes(self, uri):
        """
        Load uri contents as stream
        :param uri: URI string
        :return: content bytes
        """
        uri = os.path.join(uri)

        if os.path.isfile(uri):
            return open(uri, 'rb')
        elif os.path.isdir(uri):
            return None
