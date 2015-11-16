import os
import re
import httpretty
from tests.mock import MockResource


class GenbankMock(MockResource):
    """
    Mock Genbank resource for test purposes
    """
    def prepare_responses(self):
        """
        Prepares mock Genbank replies to requests
        """
        uris = {
            '/entrez/eutils/esearch.fcgi?.*':
                'genbank_mock_data/esearch_response.xml',
            '/entrez/eutils/efetch.fcgi?.*':
                'genbank_mock_data/mock_contents.txt'
        }
        content_types = {
            '.json': 'application/json',
            '.bam': 'application/octet-stream',
            '.bed': 'text/csv',
            '.vcf': 'text/csv',
            '.xml': 'text/xml',
            '.fasta': 'text/plain'
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
