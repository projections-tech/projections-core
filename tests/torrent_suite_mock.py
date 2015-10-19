import os
import re
import json
import httpretty
from tests.mock import MockResource

class TorrentSuiteMock(MockResource):

    def mock_auth_response(self):
        """
        Defines mock resource response to authorization request
        """
        uri = 'http://{}/'.format(self.mock_url)
        httpretty.register_uri(httpretty.GET, uri=uri, status=200)

    def prepare_responses(self):
        """
        Prepares mock Torrent Suite replies to requests
        """
        uri_dict = {
            '/rundb/api/v1/experiment\?status=run\&limit=(\d+)\&order_by=-id':
                'experiments.json',
            '/rundb/api/v1/experiment/24/':
                'experiments_24.json',
            '/rundb/api/v1/experiment/25/':
                'experiments_25.json',
            '/rundb/api/v1/plannedexperiment/31/':
                'plannedexperiment_31.json',
            '/rundb/api/v1/plannedexperiment/32/':
                'plannedexperiment_32.json',
            '/rundb/api/v1/results/8/':
                'results_8.json',
            '/rundb/api/v1/results/9/':
                'results_9.json',
            '/rundb/api/v1/pluginresult/(\d+)/':
                'plugin_result.json',
            '/rundb/api/v1/sample/(\d+)/':
                'sample_1.json',
            '/auth/output/Home/.*/plugin_out/variantCaller_out[\.\d+]*/local_parameters.json':
                'mock_vc_parameters.json',
            '/auth/output/Home/.*/IonXpress_00\d+_rawlib.bam':
                'mock_bam.bam',
            '/auth/output/Home/.*/plugin_out/variantCaller_out[\.\d+]*/.*bed':
                'mock_bed.bed',
            '/auth/output/Home/test_run/plugin_out/variantCaller_out[\.\d+]*/IonXpress_00(\d+)/.*vcf':
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

    def get_experiments(self):
        """
        Returns list of current mock experiments names
        :return: list of experiments names strings
        """
        with open(os.path.join(self.content_dir, 'experiments.json')) as ex_f:
            experiments = json.load(ex_f)
            return [exp['displayName'] for exp in experiments['objects']]

