#!/usr/bin/env python3

import json
import logging
import logging.config
import os

from projections import ProjectionDriver

logger = logging.getLogger('filesystem_projection')


class FSDriver(ProjectionDriver):
    def __init__(self, uri, driver_config, script_dir):
        self.uri = uri
        self.driver_config = driver_config
        self.script_dir = script_dir


    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """
        uri = os.path.join(self.script_dir, uri)

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
        uri = os.path.join(self.script_dir, uri)

        if os.path.isfile(uri):
            with open(uri, 'rb') as f:
                return f.read()
        elif os.path.isdir(uri):
            return json.dumps(os.listdir(uri)).encode()
