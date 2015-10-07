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
        """
        last_request = httpretty.last_request()
        return httpretty.last_request().method, last_request.path

    def mock_auth_response(self):
        """
        Reply to basic authorization call, should be overridden in resource-specific manner
        """
        pass

    def prepare_responses(self):
        """
        Prepare responses for resource specific set of URI's. Should be overridden in resource-specific manner
        """
        pass

    def __del__(self):
        httpretty.disable()
        httpretty.reset()


class TorrentSuiteMock(MockResource):

    def mock_auth_response(self):
        uri = 'http://{}/rundb/api/v1/'.format(self.mock_url)
        httpretty.register_uri(httpretty.GET, uri=uri, status=200)

    def prepare_responses(self):
        """
        Prepares mock Torrent Suite replies to requests
        """
        uri_dict = {
            '/rundb/api/v1/experiment\?status=run\&limit=1\&order_by=-id':
                'experiments_metadata_28.json',
            '/rundb/api/v1/experiment/(\d+)':
                'experiments_metadata_28.json',
            '/rundb/api/v1/plannedexperiment/(\d+)/':
                'plannedexperiment_metadata_31.json',
            '/rundb/api/v1/results/(\d+)/':
                'results_8.json',
            '/rundb/api/v1/pluginresult\?result\=(\d+)':
                'plugin_result.json',
            '/rundb/api/v1/sample/(\d+)/':
                'plugin_result.json',
            '/auth/output/Home/Run_11_hg19_v3_008/plugin_out/variantCaller_out[\.\d+]*/local_parameters.json':
                'mock_vc_parameters.json',
            '/auth/output/Home/Run_11_hg19_v3_008/IonXpress_00\d+_rawlib.bam':
                'mock_bam.bam',
            '/auth/output/Home/Run_11_hg19_v3_008/plugin_out/variantCaller_out[\.\d+]*/IAD39777_BED_4_for_TSVC.bed':
                'mock_bed.bed',
            '/auth/output/Home/Run_11_hg19_v3_008/plugin_out/variantCaller_out[\.\d+]*/IonXpress_00\d+/*.vcf':
                'mock_vcf.vcf'
        }
        content_types_dict = {
            '.json': 'application/json',
            '.bam': 'application/bam',
            '.bed': 'application/bed',
            '.vcf': 'application/vcf'
        }
        for uri, file_name in uri_dict.items():
            _, file_extension = os.path.splitext(file_name)
            if file_extension in content_types_dict:
                content_type = content_types_dict[file_extension]
            else:
                content_type = ''
            with open(os.path.join(self.content_dir, file_name), 'rb') as f:
                httpretty.register_uri(httpretty.GET,
                                       re.compile(self.mock_url + uri),
                                       body=f.read(),
                                       status=200,
                                       content_type=content_type,
                                       match_querystring=True)