__author__ = 'vsvekolkin'

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
        self.mock_url = self.mock.pretend_url

    def get_requests_to_mock(self):
        return self.mock.get_request(0)

    def mock_auth_response(self):
        pass

    def prepare_responses(self):
        self.mock.when('GET /hello').reply('Hello!'.encode(), times=FOREVER)