__author__ = 'vsvekolkin'

import os
import json
import logging
import logging.config
import urllib.request
from urllib.parse import urljoin, urlparse
from pretenders.client.http import HTTPMock
from pretenders.common.constants import FOREVER

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('resource_pretender')


class Mock_Resource(object):
    """
    Basic mock resource class to build up upon. Must be subclassed in application specific manner.
    """
    def __init__(self, content_dir):
        # Assume a running server
        # Initialise the mock client and clear all responses
        self.mock = HTTPMock('localhost', 8000)
        parse_mock_url = urlparse(self.mock.pretend_url)[1:3]
        self.mock_url = ''.join(parse_mock_url)
        logger.debug('Mock resource URL: %s', self.mock_url)
        self.content_dir = content_dir
        self.mock_auth_response()
        self.prepare_responses()
        logger.debug('Mock resource contents: %s', os.listdir(content_dir))

    def get_requests_to_mock(self):
        return self.mock.get_request(0)

    def mock_auth_response(self):
        return self.mock.reply(status=200, times=1)

    def prepare_responses(self):
        pass


class Torrent_Suite_Mock(Mock_Resource):
    def prepare_responses(self):
        rules_dict = {
            'GET /rundb/api/v1/experiment\?status\=run\&limit\=1\&order_by\=-id': 'experiments_metadata.json',
            'GET /rundb/api/v1/plannedexperiment/31/': 'plannedexperiment_metadata.json',
            'GET /rundb/api/v1/results/8/': 'results.json',
            'GET /rundb/api/v1/pluginresult\?result\=8': 'plugin_result.json'
        }
        for rule, json_file_name in rules_dict.items():
            with open(os.path.join(self.content_dir, json_file_name), 'rb') as f:
                self.mock.when(rule).reply(body=f.read(),
                                           headers={'Content-Type': 'application/json'},
                                           times=FOREVER)
        for i in range(6):
            rule = 'GET /auth/output/Home/Run_11_hg19_v3_008/IonXpress_00{0}_rawlib.bam'.format(i)
            self.mock.when(rule).reply(body=b'Mock BAM file here.', times=1)
            for vc_dir in ['variantCaller_out', 'variantCaller_out.49', 'variantCaller_out.50']:
                for variant_file_name in ['TSVC_variants.vcf', 'all.merged.vcf', 'indel_assembly.vcf',
                                          'indel_variants.vcf', 'small_variants.left.vcf',
                                          'small_variants.vcf', 'small_variants_filtered.vcf',
                                          'small_variants.sorted.vcf', 'SNP_variants.vcf']:
                    rule = 'GET /auth/output/Home/Run_11_hg19_v3_008/plugin_out/{0}/IonXpress_00{1}/{2}'.format(vc_dir,
                                                                                                                i, variant_file_name)
                    self.mock.when(rule).reply(body=b'Mock VCF file here.', times=1)


if __name__ == '__main__':
    mocker = Torrent_Suite_Mock('mock_resource')
    url_p = 'http://{}/rundb/api/v1/experiment?status=run&limit=1&order_by=-id'.format(mocker.mock_url)
    logger.debug(url_p)
    with urllib.request.urlopen(url_p) as f:
        experiments = json.loads(f.readall().decode('utf-8'))

