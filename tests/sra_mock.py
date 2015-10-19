import os
import re
import httpretty
from tests.mock import MockResource


class SRAMock(MockResource):
    def prepare_responses(self):
        """
        Prepares mock SRA replies to requests
        """
        uri_dict = {
            '/entrez/eutils/esearch.fcgi.*':
                'esearch_query.xml',
            '/entrez/eutils/efetch.fcgi.*':
                'efetch_query.xml'
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
