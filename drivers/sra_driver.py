#!/usr/bin/env python3

#    Copyright 2016  Anton Bragin, Victor Svekolkin
#
#    This file is part of Projections.
#
#    Projections is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Projections is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Projections.  If not, see <http://www.gnu.org/licenses/>.

import json
import logging
import logging.config
import re
import subprocess

import xmltodict
from Bio import Entrez

from projections import ProjectionDriver

logger = logging.getLogger('sra_driver')


class SRADriver(ProjectionDriver):
    def __init__(self, uri, driver_config_path):
        """
        Initialize driver telling NCBI using user email and name of program
        :param driver_config: driver configuration dictionary
        :return:
        """
        self.driver_configuration = self.read_config(driver_config_path)

        Entrez.email = self.driver_configuration['email']
        Entrez.tool = 'sra_projection_manager'

        self.sampdump_path = self.driver_configuration['samdump_path']

    def get_uri_contents_as_dict(self, uri):
        """
        Loads content from SRA using Biopython in driver cache
        :param uri: str containing uri
        :return: dict of query contents
        """

        sam_uri_regex = '(SRR|SRX|ERX|DRX|DRR|ERR)\d+'
        if re.match(sam_uri_regex, uri):
            return {}

        if uri.startswith('search_query:'):
            query = uri.split(':')

            if len(query) < 3:
                query.append(1)
            esearch_handle = Entrez.esearch(db='sra', term=query[1], retmax=query[2])

            search_result = Entrez.read(esearch_handle)
            search_result['IdList'] = ['sra_id:' + el for el in search_result['IdList']]

            return dict(search_result)
        elif uri.startswith('sra_id:'):
            query = uri.split(':')
            fetch_handler = Entrez.efetch(db='sra', id=query[1])

            experiment = json.loads(json.dumps(xmltodict.parse(fetch_handler.read())))

            # Experiment resource contains data about runs, which contain id`s of SAM files in database
            # Run set is most times dict, but sometimes list, treating dict as list to resolve inconsistency
            run_set = experiment['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN']
            if not isinstance(run_set, list):
                run_set = [run_set]
            # Setting corrected run set field for experiment
            experiment['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN'] = run_set
            return dict(experiment)

    def get_uri_contents_as_bytes(self, uri):
        """
        Returns stream of uri contents
        :param uri: uri of resource
        :return: stream of resource contents
        """
        logger.debug('Loading query: %s', uri)
        # Regex matches id`s in SRA database
        sam_uri_regex = '(SRR|SRX|ERX|DRX|DRR|ERR)\d+'
        if re.match(sam_uri_regex, uri):
            # Using subprocess.check to run sam-dump, which returns stream after loading of sam file,
            # this approach may be slow, other approaches lock script execution, need to reconsider
            return SamDump(self.sampdump_path, uri)
        else:
            return (el for el in json.dumps(self.get_uri_contents_as_dict(uri)).encode())


class SamDump:
    """
    This class works as context manager for sam-damp tool.
    """

    def __init__(self, sampdump_path, uri):
        self.sampdump_path = sampdump_path
        self.uri = uri
        self.sam_dump = None

    def __enter__(self):
        self.sam_dump = subprocess.Popen([self.sampdump_path, self.uri], stdout=subprocess.PIPE)
        return self.sam_dump.stdout

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self.sam_dump.terminate()
        self.sam_dump.terminate()
        self.sam_dump.communicate()
