__author__ = 'vsvekolkin'

import urllib.request
import json
import logging
from pretenders.client.http import HTTPMock
from pretenders.common.constants import FOREVER

class Mock_Resource:
    """
    Basic mock resource class to build up upon. Must be subclassed in application specific manner.
    """
    def __init__(self):
        # Assume a running server
        # Initialise the mock client and clear all responses
        self.mock = HTTPMock('localhost', 8000)
        self.mock_url = self.mock.pretend_url+'/'
        print(self.mock_url)
        self.prepare_responses()

    def get_requests_to_mock(self):
        return self.mock.get_request(0)

    def mock_auth_response(self):
        pass

    def prepare_responses(self):
        self.mock.when('GET /experiment').reply('{100:101}'.encode(), status=200, times=FOREVER)


if __name__ == '__main__':
    res_mock = Mock_Resource()
    print(urllib.parse.urljoin(res_mock.mock_url, 'experiment'))
    with urllib.request.urlopen(urllib.parse.urljoin(res_mock.mock_url, 'experiment')) as f:
        print(f.readall())