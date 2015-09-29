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
        self.mock_auth_response()
        self.prepare_responses()

    def get_requests_to_mock(self):
        return self.mock.get_request(0)

    def mock_auth_response(self):
        return self.mock.when(rule='GET /',
                              headers={'User-Agent': 'Python-urllib/3.4',
                                       'Accept-Encoding': 'identity',
                                       'Content-Type': 'text/plain',
                                       'Connection': 'close',
                                       'Host': 'localhost:8000',
                                       'Content-Length': ''}).reply(status=200)

    def prepare_responses(self):
        pass


if __name__ == '__main__':
    res_mock = Mock_Resource()
    password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    user = 'User'
    password = 'Test'
    password_manager.add_password(None, res_mock.mock_url, user, password)
    handler = urllib.request.HTTPBasicAuthHandler(password_manager)
    # create "opener" (OpenerDirector instance)
    opener = urllib.request.build_opener(handler)
    # use the opener to fetch a URL
    opener.open(res_mock.mock_url)
    # Install the opener.
    # Now all calls to urllib.request.urlopen use our opener.
    urllib.request.install_opener(opener)
