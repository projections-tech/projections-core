__author__ = 'vsvekolkin'

import os
import re
import json
import logging
import logging.config
import httpretty

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('mock')


class MockResource(object):
    """
    Mock resource
    """
    def __init__(self, configuration_path):
        """
        Initializes mock resource with provided JSON configuration file
        :param configuration_path path to configuration JSON file string
        """
        logger.debug('Initializing mock resource.')
        httpretty.enable()
        # Disable network connection for mock server, this option enables more transparent error messages
        httpretty.HTTPretty.allow_net_connect = False
        # Loading configuration of mock from JSON file
        with open(configuration_path) as m_c:
            self.mock_configuration = json.load(m_c)
        # Directory with mock contents files
        self.content_dir = self.mock_configuration['mock_data_path']
        # Url which will be mocked
        self.mock_url = self.mock_configuration['mock_resource_uri']
        logger.debug('Mock resource URL: %s', self.mock_url)
        logger.debug('HTTPpretty status: %s', httpretty.is_enabled())
        logger.debug('Mock resource contents: %s', os.listdir(self.content_dir))

        # Setting up mock authorization response
        self.mock_auth_response()
        # Registering URI`s which will be mocked
        self.prepare_responses(self.mock_configuration['mock_responses'])

    def get_last_request_to_mock(self):
        """
        Returns last request to current mock
        :return: tuple. first element is response method, second is path
        """
        return httpretty.last_request().method, httpretty.last_request().path

    def mock_auth_response(self):
        """
        Response to basic auth call
        """
        uri = 'http://{}/'.format(self.mock_url)
        httpretty.register_uri(httpretty.GET, uri=uri, status=200)

    def prepare_responses(self, uris):
        """
        Prepares responses for resource specific set of URI's, provided in config.
        :param uris dict of URI content file path pairs
        """
        content_types = {
            '.json': 'application/json',
            '.bam': 'application/octet-stream',
            '.bed': 'text/csv',
            '.vcf': 'text/csv',
            '.xml': 'text/xml'
        }

        for uri, file_name in uris.items():
            _, file_extension = os.path.splitext(file_name)
            if file_extension in content_types:
                content_type = content_types[file_extension]
            else:
                content_type = ''
            with open(os.path.join(self.content_dir, file_name), 'rb') as f:
                httpretty.register_uri(httpretty.GET,
                                       re.compile(self.mock_url + uri),
                                       body=f.read(),
                                       status=200,
                                       content_type=content_type,
                                       match_querystring=True)

    def __del__(self):
        """
        Destructor of mock resource object
        """
        httpretty.disable()
        httpretty.reset()
