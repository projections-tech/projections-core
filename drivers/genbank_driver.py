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

from Bio import Entrez

from projections import ProjectionDriver

logger = logging.getLogger('genbank_driver')


class GenbankDriver(ProjectionDriver):
    def __init__(self, uri, driver_config_path):
        self.driver_configuration = self.read_config(driver_config_path)

        Entrez.email = self.driver_configuration['email']
        Entrez.tool = 'genbank_projection_manager'
        self.driver_cache = {}

    def get_uri_contents_as_dict(self, uri):
        """
        Loads content from SRA using Biopython, with queries in format: "query_type:query"
        :param uri: str containing query to SRA
        :return: dict of query contents
        """
        if uri.startswith('gb:') or uri.startswith('fasta:'):
            return {}

        if uri not in self.driver_cache and uri.startswith('search_query:'):
            query = uri.split(':')
            logger.debug('Current query: %s', query)
            # Info about Biopython`s eutils: http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter:entrez
            # Query looks as: 'query:Test_species'

            # Returns esearch response dict.
            esearch_handle = Entrez.esearch(db='nucleotide', term=query[1], retmax=query[2])
            search_result = dict(Entrez.read(esearch_handle))
            search_result['resource_uri'] = uri
            search_result['fasta_files'] = []
            search_result['gb_files'] = []

            for nuc_id in search_result['IdList']:
                fasta_id = 'fasta:{}'.format(nuc_id)
                gb_id = 'gb:{}'.format(nuc_id)
                query_id = 'query:'.format(nuc_id)

                search_result['fasta_files'].append(fasta_id)
                search_result['gb_files'].append(gb_id)

                self.driver_cache[fasta_id] = {'resource_uri': fasta_id}
                self.driver_cache[gb_id] = {'resource_uri': gb_id}
                self.driver_cache[query_id] = {'resource_uri': query_id}
            self.driver_cache[uri] = search_result
        return self.driver_cache[uri]

    def get_uri_contents_as_bytes(self, query):
        """
        Open URI and return it`s content as bytes
        :param query: query to driver
        :return: bytes massive
        """
        logger.debug('Loading query: %s', query)
        query = query.split(':')
        query_type = query[0]
        if query_type == 'search_query':
            return Entrez.esearch(db='nuccore', term=query[1], retmax=query[2]).read().encode()
        elif query_type == 'gb':
            # Returns query gb file bytes
            gb = Entrez.efetch(db='nuccore', id=query[1], rettype='gb', retmode='text')
            return gb.read().encode()
        elif query_type == 'fasta':
            # Returns query fasta file bytes
            fasta = Entrez.efetch(db='nuccore', id=query[1], rettype='fasta', retmode='text')
            return fasta.read().encode()
