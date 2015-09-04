__author__ = 'abragin'

import logging
import json
import urllib.request
from urllib.parse import urljoin
from projections import ProjectionManager


logger = logging.getLogger('iontorrent_projection')

class IonTorrentProjection(ProjectionManager):

    def __init__(self, host, user, password):
        self.api_url = 'http://{}/rundb/api/v1/'.format(host)
        self.authenticate(user, password)

    def authenticate(self, user, password):
        password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        password_manager.add_password(None, self.api_url, user, password)

        handler = urllib.request.HTTPBasicAuthHandler(password_manager)

        # create "opener" (OpenerDirector instance)
        opener = urllib.request.build_opener(handler)

        # use the opener to fetch a URL
        opener.open(self.api_url)

        # Install the opener.
        # Now all calls to urllib.request.urlopen use our opener.
        urllib.request.install_opener(opener)


    def create_projections(self):
        with urllib.request.urlopen(urljoin(self.api_url, 'experiment/')) as f:
            experiments = json.loads(f.readall().decode('utf-8'))
        logger.info('Got experiments data: %s', experiments)






