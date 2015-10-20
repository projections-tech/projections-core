import os
import re
import httpretty
from tests.mock import MockResource


class SRAMock(MockResource):
    """
    Mock SRA resource for test purposes
    """
    def prepare_responses(self):
        """
        Prepares mock SRA replies to requests
        """
        uris = {
            '/entrez/eutils/esearch.fcgi.*':
                'sra_mock_data/esearch_query.xml',
            '/entrez/eutils/efetch.fcgi.*':
                'sra_mock_data/efetch_query.xml'
        }
        content_types = {
            '.json': 'application/json',
            '.bam': 'application/octet-stream',
            '.bed': 'text/csv',
            '.vcf': 'text/csv',
            '.xml': 'text/xml; charset=UTF-8'
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
