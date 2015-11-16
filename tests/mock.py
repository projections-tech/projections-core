__author__ = 'vsvekolkin'

import os
import re
import logging
import logging.config
import httpretty

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('mock')


class MockResource(object):
    """
    Basic mock resource class to build up upon. Must be subclassed in application specific manner.
    """
    def __init__(self, resource_url, content_dir):
        """
        Initializes mock resource
        """
        logger.debug('Initializing mock resource.')
        httpretty.enable()
        # Disable network connection for mock server, this option enables more transparent error messages
        httpretty.HTTPretty.allow_net_connect = False
        self.content_dir = content_dir
        self.mock_url = resource_url
        logger.debug('Mock resource URL: %s', self.mock_url)
        logger.debug('HTTPpretty status: %s', httpretty.is_enabled())
        logger.debug('Mock resource contents: %s', os.listdir(content_dir))
        # Setting up mock authorization response
        self.mock_auth_response()
        self.prepare_responses()

    def get_last_request_to_mock(self):
        """
        Returns last request to current mock
        :return: tuple. first element is response method, second is path
        """
        return httpretty.last_request().method, httpretty.last_request().path

    def mock_auth_response(self):
        """
        Reply to basic authorization call, should be overridden in resource-specific manner
        """
        return NotImplementedError('Authorization is not implemented!')

    def prepare_responses(self):
        """
        Prepare responses for resource specific set of URI's. Should be overridden in resource-specific manner
        """
        return NotImplementedError('Responses are not implemented in base class, this method must be implemented in subclass in application-specific manner.')

    def deactivate(self):
        httpretty.disable()
        httpretty.reset()
