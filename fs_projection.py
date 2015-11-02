#!/usr/bin/env python3

import logging
import logging.config
import io
import json
import os
import sys
import time

from projections import Projection, ProjectionDriver, ProjectionTree, Projector, PrototypeDeserializer
from filesystem import ProjectionFilesystem
from fuse import FUSE

logger = logging.getLogger('filesystem_projection')


class FSDriver(ProjectionDriver):
    def __init__(self, host_uri):
        """
        Initialize driver which will be used to interact with host.
        """
        self.host_uri = host_uri

    def get_uri_contents_as_dict(self, uri):
        """
        Opens URI and returns dict of its contents
        :param uri: URI string
        :return: dict of URI contents
        """
        if os.path.isdir(uri):
            return json.dumps(os.listdir(uri))
        else:
            return json.dumps({'name':os.path.basename(uri), 'size': os.path.getsize(uri),
                               'extension': os.path.splitext(uri)[1]})

    def get_uri_contents_as_stream(self, uri):
        """
        Load uri contents
        :param uri: URI string
        :return: content bytes
        """
        if os.path.isfile(uri):
            return open(uri, 'rb')
        elif os.path.isdir(uri):
            return json.dumps(os.listdir(uri)).encode()


if __name__ == '__main__':
    driver = FSDriver('tests')
    print(driver.get_uri_contents_as_stream('tests'))
