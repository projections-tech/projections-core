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

import logging
import logging.config
import os

import requests
from Bio import Entrez

from projections import ProjectionDriver

logger = logging.getLogger('genbank_driver')


class GenbankDriver(ProjectionDriver):
    def __init__(self, uri, driver_config_path):
        self.driver_configuration = self.read_config(driver_config_path)

        Entrez.email = self.driver_configuration['email']
        Entrez.tool = 'genbank_projection_manager'

    def get_uri_contents_as_dict(self, uri):
        """
        Loads content from SRA using Biopython, with queries in format: "query_type:query"
        :param uri: str containing query to SRA
        :return: dict of query contents
        """

        uri, _ = os.path.splitext(uri)

        if uri.startswith('search_query:'):
            query = uri.split(':')
            logger.debug('Current query: %s', query)
            if len(query) < 3:
                query.append(1)
            # Info about Biopython`s eutils: http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter:entrez
            # Returns esearch response dict.
            esearch_handle = Entrez.esearch(db='nucleotide', term=query[1], retmax=query[2])
            return dict(Entrez.read(esearch_handle))
        else:
            esearch_handle = Entrez.esearch(db='nucleotide', term=uri, retmax=1)
            return dict(Entrez.read(esearch_handle))

    def get_uri_contents_as_bytes(self, query):
        """
        Open URI and return it`s content as bytes
        :param query: query to driver
        :return: bytes massive
        """
        query, query_type = os.path.splitext(query)

        logger.debug('Loading query: %s', query)
        if query_type == '.json':
            return (el for el in Entrez.esearch(db='nuccore', term=query[1], retmax=query[2]).read().encode())
        elif query_type == '.gb':
            # Returns query gb file bytes iterator
            return self.load_data_from_genbank(query, 'gb')
        elif query_type == '.fasta':
            # Returns query fasta file bytes iterator
            return self.load_data_from_genbank(query, 'fasta')

    def load_data_from_genbank(self, req_id, ret_mode):
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' \
              '?db=nucleotide&id={acc}&rettype={ret_type}&retmode=text"'.format(acc=req_id,
                                                                                ret_type=ret_mode)
        r = requests.get(url, stream=True)
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:
                yield chunk
