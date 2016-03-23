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
    def __init__(self, uri, driver_config_path, script_dir):
        """
        Initialize driver telling NCBI using user email and name of program
        :param driver_config: driver configuration dictionary
        :return:
        """
        self.driver_configuration = self.read_config(script_dir, driver_config_path)

        Entrez.email = self.driver_configuration['email']
        Entrez.tool = 'sra_projection_manager'
        self.driver_cache = {}

    def get_uri_contents_as_dict(self, uri):
        """
        Loads content from SRA using Biopython in driver cache
        :param uri: str containing uri
        :return: dict of query contents
        """

        sam_uri_regex = '(SRR|SRX|ERX|DRX|DRR|ERR)\d+'
        if re.match(sam_uri_regex, uri):
            return {}

        # Info about Biopython`s eutils: http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter:entrez
        # Query looks as: 'query:Test_species'
        logger.debug('Current query: %s', uri)
        uri_parts = uri.split(':')
        if uri not in self.driver_cache and uri_parts[0] == 'search_query':
            # Returns esearch response dict for SRA database.
            logger.debug('Cureent uri parts: %s', uri_parts)
            esearch_handle = Entrez.esearch(db='sra', term=uri_parts[1], retmax=uri_parts[2])
            search_result = Entrez.read(esearch_handle)
            # Id`s of experiment resources
            for sra_id in search_result['IdList']:
                # Driver accesses experiments using id`s like sra_id:1214564
                res_uri = 'sra_id:{0}'.format(sra_id)
                # Handler to fetch SRA database by experiment id
                fetch_handler = Entrez.efetch(db='sra', id=sra_id)

                # Converts nested ordered dicts to default dicts required by ObjectPath
                experiment = json.loads(json.dumps(xmltodict.parse(fetch_handler.read())))
                # Setting resource uri for experiment on driver
                experiment['resource_uri'] = res_uri

                # Experiment resource contains data about runs, which contain id`s of SAM files in database
                # Run set is most times dict, but sometimes list, treating dict as list to resolve inconsistency
                run_set = experiment['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN']
                if not isinstance(run_set, list):
                    run_set = [run_set]
                # Setting corrected run set field for experiment
                experiment['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN'] = run_set
                # Adding experiment resource in driver cache by uri
                self.driver_cache[res_uri] = experiment

                # Adding runs to driver cache by their accession
                for run in run_set:
                    self.driver_cache[run['@accession']] = run

            # Adding sra_id prefix by which experiment resources will be accessed using driver cache
            search_result['IdList'] = [''.join(['sra_id:', id]) for id in search_result['IdList']]
            # Adding resource uri to search results meta
            search_result['resource_uri'] = uri
            # Converting search result into proper dict, not biopython subclass, because ObjectPath can handle only
            # standard dicts
            search_result = dict(search_result)
            # Adding results of search to driver cache
            self.driver_cache[uri] = search_result
            return search_result
        elif uri_parts[0] == 'sra_id':
            # Driver accesses experiments using id`s like sra_id:1214564
            res_uri = 'sra_id:{0}'.format(uri_parts[1])
            # Handler to fetch SRA database by experiment id
            fetch_handler = Entrez.efetch(db='sra', id=uri_parts[1])

            # Converts nested ordered dicts to default dicts required by ObjectPath
            experiment = json.loads(json.dumps(xmltodict.parse(fetch_handler.read())))
            # Setting resource uri for experiment on driver
            experiment['resource_uri'] = res_uri

            # Experiment resource contains data about runs, which contain id`s of SAM files in database
            # Run set is most times dict, but sometimes list, treating dict as list to resolve inconsistency
            run_set = experiment['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN']
            if not isinstance(run_set, list):
                run_set = [run_set]
            # Setting corrected run set field for experiment
            experiment['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN'] = run_set
            # Adding experiment resource in driver cache by uri
            self.driver_cache[res_uri] = experiment

            # Adding runs to driver cache by their accession
            for run in run_set:
                self.driver_cache[run['@accession']] = run
            return experiment
        else:
            return self.driver_cache[uri]

    def load_uri_contents_as_bytes(self, uri):
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
            return subprocess.check_output(['./sratoolkit.2.5.4-1-ubuntu64/bin/sam-dump', uri])
        else:
            return json.dumps(self.driver_cache[uri]).encode()