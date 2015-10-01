__author__ = 'vsvekolkin'

import os
import re
import logging
import logging.config
import httpretty
import urllib

logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('resource_pretender')


class MockResource(object):
    """
    Basic mock resource class to build up upon. Must be subclassed in application specific manner.
    """
    def __init__(self, content_dir):
        """
        Initializes mock resource
        """
        logger.debug('Initializing mock resource.')
        httpretty.enable()
        # Disable network connection for mock server, this option enables more transparent error messages
        httpretty.HTTPretty.allow_net_connect = False
        self.content_dir = content_dir
        self.mock_url = 'torrentsuitemock.com'
        logger.debug('Mock resource URL: %s', self.mock_url)
        logger.debug('HTTPpretty status: %s', httpretty.is_enabled())
        logger.debug('Mock resource contents: %s', os.listdir(content_dir))
        # Setting up mock authorization response
        self.mock_auth_response()
        self.prepare_responses()

    def get_requests_to_mock(self):
        """
        Returns list of all requests to current mock
        """
        return httpretty.last_request()

    def mock_auth_response(self):
        """
        Reply to basic authorization call
        """
        uri = 'http://{}/rundb/api/v1/'.format(self.mock_url)
        httpretty.register_uri(httpretty.GET, uri=uri, status=200)

    def prepare_responses(self):
        """
        This method should be overridden in resource-specific manner
        """
        pass

    def deactivate_mock(self):
        httpretty.disable()
        httpretty.reset()


class TorrentSuiteMock(MockResource):
    def prepare_responses(self):
        """
        Prepares mock Torrent Suite replies to requests
        """
        rules_dict = {
            '/rundb/api/v1/experiment\?status\=run\&limit\=1\&order_by\=-id': 'experiments_metadata.json',
            '/rundb/api/v1/experiment/\d+': 'experiments_metadata.json',
            '/rundb/api/v1/plannedexperiment/\d+/': 'plannedexperiment_metadata.json',
            '/rundb/api/v1/results/\d+/': 'results.json',
            '/rundb/api/v1/pluginresult\?result\=\d+': 'plugin_result.json',
            '/rundb/api/v1/sample/\d+/': 'plugin_result.json'
        }
        for rule, json_file_name in rules_dict.items():
            with open(os.path.join(self.content_dir, json_file_name), 'rb') as f:
                httpretty.register_uri(httpretty.GET,
                                       re.compile(self.mock_url + rule),
                                       body=f.read(),
                                       status=200,
                                       content_type='application/json', match_querystring=True)

        rule = '/auth/output/Home/Run_11_hg19_v3_008/IonXpress_00\d+_rawlib.bam'
        httpretty.register_uri(httpretty.GET, body=b'Mock BAM file here.',
                               uri=re.compile(self.mock_url + rule), status=200)

        rule = '/auth/output/Home/Run_11_hg19_v3_008/plugin_out/variantCaller_out[\.\d+]*/IAD39777_BED_4_for_TSVC.bed'
        httpretty.register_uri(httpretty.GET, body=b'Mock BED file here.',
                       uri=re.compile(self.mock_url + rule), status=200)

        rule = '/auth/output/Home/Run_11_hg19_v3_008/plugin_out/variantCaller_out[\.\d+]*/local_parameters.json'
        httpretty.register_uri(httpretty.GET, body=b'Mock VC settings file here.',
                               uri=re.compile(self.mock_url + rule), status=200)

        for variant_file_name in ['TSVC_variants.vcf', 'all.merged.vcf', 'indel_assembly.vcf',
                                  'indel_variants.vcf', 'small_variants.left.vcf',
                                  'small_variants.vcf', 'small_variants_filtered.vcf',
                                  'small_variants.sorted.vcf', 'SNP_variants.vcf']:
            rule = '/auth/output/Home/Run_11_hg19_v3_008/plugin_out/variantCaller_out[\.\d+]*/IonXpress_00\d+/{0}'.format(re.escape(variant_file_name))
            httpretty.register_uri(httpretty.GET, body=b'Mock VCF file here.',
                                   uri=re.compile(self.mock_url + rule), status=200)

class SRAMock(MockResource):
    def prepare_responses(self):
        """
        Prepares mock SRA replies to requests
        """
        pass
