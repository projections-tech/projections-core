__author__ = 'vsvekolkin'

import os
import re
import json
import logging
import logging.config
import urllib.request
from urllib.parse import urljoin, urlparse
from pretenders.client.http import HTTPMock
from pretenders.common.constants import FOREVER

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('resource_pretender')


class MockResource(object):
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


class TorrentSuiteMock(MockResource):
    def prepare_responses(self):
        rules_dict = {
            'GET /rundb/api/v1/experiment\?status\=run\&limit\=1\&order_by\=-id': 'experiments_metadata.json',
            'GET /rundb/api/v1/experiment/\d+': 'experiments_metadata.json',
            'GET /rundb/api/v1/plannedexperiment/\d+/': 'plannedexperiment_metadata.json',
            'GET /rundb/api/v1/results/\d+/': 'results.json',
            'GET /rundb/api/v1/pluginresult\?result\=\d+': 'plugin_result.json',
            'GET /rundb/api/v1/sample/\d+/': 'plugin_result.json'
        }
        for rule, json_file_name in rules_dict.items():
            with open(os.path.join(self.content_dir, json_file_name), 'rb') as f:
                self.mock.when(rule).reply(body=f.read(),
                                           headers={'Content-Type': 'application/json'},
                                           times=FOREVER)

        rule = 'GET /auth/output/Home/Run_11_hg19_v3_008/IonXpress_00\d+_rawlib.bam'
        self.mock.when(rule).reply(body=b'Mock BAM file here.', times=FOREVER, status=200)

        rule = 'GET /auth/output/Home/Run_11_hg19_v3_008/plugin_out/variantCaller_out[\.\d+]*/IAD39777_BED_4_for_TSVC.bed'
        self.mock.when(rule).reply(body=b'Mock BED file here.', times=FOREVER, status=200)

        rule = 'GET /auth/output/Home/Run_11_hg19_v3_008/plugin_out/variantCaller_out[\.\d+]*/local_parameters.json'
        self.mock.when(rule).reply(body=b'Mock VC settings file here.', times=FOREVER, status=200)

        for variant_file_name in ['TSVC_variants.vcf', 'all.merged.vcf', 'indel_assembly.vcf',
                                  'indel_variants.vcf', 'small_variants.left.vcf',
                                  'small_variants.vcf', 'small_variants_filtered.vcf',
                                  'small_variants.sorted.vcf', 'SNP_variants.vcf']:
            rule = 'GET /auth/output/Home/Run_11_hg19_v3_008/plugin_out/variantCaller_out[\.\d+]*/IonXpress_00\d+/{0}'.format(re.escape(variant_file_name))
            self.mock.when(rule).reply(body=b'Mock VCF file here.', status=200, times=FOREVER)


if __name__ == '__main__':
    mocker = TorrentSuiteMock('mock_resource')
    url_p = 'http://{}/rundb/api/v1/experiment?status=run&limit=1&order_by=-id'.format(mocker.mock_url)
    logger.debug(url_p)
    with urllib.request.urlopen(url_p) as f:
        experiments = json.loads(f.readall().decode('utf-8'))

