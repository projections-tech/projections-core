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
        uris = {
            '/rundb/api/v1/experiment\?status=run\&limit=(\d+)\&order_by=-id':
                'experiments.json',
            '/rundb/api/v1/experiment/1/':
                'experiments_1.json',
            '/rundb/api/v1/experiment/2/':
                'experiments_2.json',
            '/rundb/api/v1/plannedexperiment/1/':
                'plannedexperiment_1.json',
            '/rundb/api/v1/plannedexperiment/2/':
                'plannedexperiment_2.json',
            '/rundb/api/v1/results/1/':
                'results_1.json',
            '/rundb/api/v1/results/2/':
                'results_2.json',
            '/rundb/api/v1/results/3/':
                'results_3.json',
            '/rundb/api/v1/results/4/':
                'results_4.json',
            '/rundb/api/v1/pluginresult/(\d+)/':
                'plugin_result.json',
            '/rundb/api/v1/sample/1/':
                'sample_1.json',
            '/rundb/api/v1/sample/2/':
                'sample_2.json',
            '/rundb/api/v1/sample/3/':
                'sample_3.json',
            '/rundb/api/v1/sample/4/':
                'sample_4.json',
            '/auth/output/Home/.*/plugin_out/variantCaller_out[\.\d+]*/local_parameters.json':
                'mock_vc_parameters.json',
            '/auth/output/Home/.*/IonXpress_00\d+_rawlib.bam':
                'mock_bam.bam',
            '/auth/output/Home/.*/plugin_out/variantCaller_out[\.\d+]*/.*bed':
                'mock_bed.bed',
            '/auth/output/Home/test_run_(\d+)/plugin_out/variantCaller_out[\.\d+]*/IonXpress_00(\d+)/.*vcf':
                'mock_vcf.vcf'
        }
        content_types = {
            '.json': 'application/json',
            '.bam': 'application/octet-stream',
            '.bed': 'text/csv',
            '.vcf': 'text/csv',
            '.xml': 'text/xml'
        }
        for uri, file_name in uris.items():
            _, file_extension = os.path.splitext(file_name)
            if file_extension in content_types:
                content_type = content_types[file_extension]
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

