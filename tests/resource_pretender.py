__author__ = 'vsvekolkin'

import os
import io
import re
import json
import logging
import logging.config
import urllib.request
from urllib.parse import urljoin
import requests
from pretenders.client.http import HTTPMock
from pretenders.common.constants import FOREVER

logging.config.fileConfig('../logging.cfg')
logger = logging.getLogger('resource_pretender')

class Mock_Resource(object):
    """
    Basic mock resource class to build up upon. Must be subclassed in application specific manner.
    """
    def __init__(self, content_dir):
        # Assume a running server
        # Initialise the mock client and clear all responses
        self.mock = HTTPMock('localhost', 8000)
        parse_mock_url = urllib.parse.urlparse(self.mock.pretend_url)[1:3]
        self.mock_url = ''.join(parse_mock_url)
        logger.debug('Mock resource URL: %s',self.mock_url)
        self.content_dir = content_dir
        self.mock_auth_response()
        self.prepare_responses()
        logger.debug('Mock resource contents: %s', os.listdir(content_dir))

    def get_requests_to_mock(self):
        return self.mock.get_request(0)

    def mock_auth_response(self):
        return self.mock.when(rule='GET /',
                              headers={'User-Agent': 'Python-urllib/3.4',
                                       'Accept-Encoding': 'identity',
                                       'Content-Type': 'text/plain',
                                       'Connection': 'close',
                                       'Host': 'localhost:8000',
                                       'Content-Length': ''}).reply(status=200, times=FOREVER)

    def prepare_responses(self):
        pass


class Torrent_Suite_Mock(Mock_Resource):
    def prepare_responses(self):
        rule = re.escape('GET /rundb/api/v1/experiment?status=run&limit=1&order_by=-id')
        logger.debug('Experiments rule: %s',rule)
        with open(os.path.join(self.content_dir, 'experiments_metadata.json'), 'r') as e_m:
            experiment_json = e_m.read()
            self.mock.when(rule).reply(body=experiment_json.encode(),
                                       headers={'Content-Type': 'application/json'}, status=200, times=3)



if __name__=='__main__':
    mocker = Torrent_Suite_Mock('mock_resource')
    url_p = 'http://'+mocker.mock_url+'/rundb/api/v1/experiment?status=run&limit=1&order_by=-id'
    logger.debug(url_p)
    print(requests.get(url_p).json())
    with urllib.request.urlopen(url_p) as f:
        logger.debug(f.readall())
        experiments = json.loads(f.readall().decode('utf-8'))