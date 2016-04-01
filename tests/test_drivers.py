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
import os
import re
from unittest import TestCase, skip

import boto3
from moto import mock_s3

from drivers.aws_s3_driver import S3Driver
from drivers.fs_driver import FSDriver
from drivers.genbank_driver import GenbankDriver
from drivers.iontorrent_driver import TorrentSuiteDriver
from drivers.sra_driver import SRADriver
from tests.mock import MockResource

logging.config.fileConfig('logging.cfg')


class TestFsDriver(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('fs_driver_test')

        cls.driver = FSDriver(driver_config_path='tests/driver_configurations/fs_config.yaml',
                              uri='')

    def test_get_uri_contents_as_dict(self):
        """
        Tests driver`s get_uri_contents_as_dict method on dirs and files
        """

        # Checking driver response for directories uri`s
        objects_paths = [os.path.join('tests/test_dir', path)
                         for path in os.listdir('tests/test_dir')]

        for path in objects_paths:
            object_name = os.path.basename(path)
            resource_uri = os.path.abspath(path)
            object_size = os.path.getsize(resource_uri)

            driver_response = self.driver.get_uri_contents_as_dict(path)

            if os.path.isdir(path):
                reference_type = 'dir'
                reference_extension = None
                num_children = len(os.listdir(resource_uri))
                # Checking children field type and children number for dir
                self.assertIsInstance(driver_response['children'], list,
                                      msg='Checking if response field "children" is correct.')
                self.assertEqual(num_children, len(driver_response['children']),
                                 msg='Checking if number of dirs children is correct.')
            else:
                reference_type = 'file'
                reference_extension = os.path.splitext(object_name)[1]
                # Checking if field "children" is in response, for files it must not be present
                self.assertNotIn('children', driver_response, msg='Checking if response field "children" is exists.')

            self.assertIsInstance(driver_response, dict, msg='Checking response object type.')
            self.assertEqual(reference_type, driver_response['type'], msg='Checking if response field type is correct.')
            self.assertEqual(object_size, driver_response['size'], msg='Checking if dir size is correct.')
            self.assertEqual(object_name, driver_response['name'], msg='Checking if dir name is correct.')
            self.assertEqual(resource_uri, driver_response['resource_uri'], msg='Checking if resource uri is correct.')
            self.assertEqual(driver_response['extension'], reference_extension, msg='Checking if extension is correct.')

    def test_get_uri_contents_as_bytes(self):
        """
        Tests driver get_uri_contents_as_bytes method on dirs and files
        """
        # Checking driver response for directories uri`s
        dirs_paths = [os.path.join('tests/test_dir', path)
                      for path in os.listdir('tests/test_dir')
                      if os.path.isdir(os.path.join('tests/test_dir', path))]

        for path in dirs_paths:
            with self.driver.get_uri_contents_as_bytes(path) as driver_response:
                response = driver_response

            self.assertIsNone(response, msg='Checking dir contents.')

        # Checking driver responses for files uri`s
        files_paths = [os.path.join('tests/test_dir', path)
                       for path in os.listdir('tests/test_dir')
                       if os.path.isfile(os.path.join('tests/test_dir', path))]

        for path in files_paths:
            # Reading reference file from test_dir
            with open(path, 'rb') as tested_file:
                reference_file_response_contents = tested_file.read()

            driver_response = self.driver.get_uri_contents_as_bytes(path)
            from drivers.fs_driver import FSLoadData
            self.assertIsInstance(driver_response, FSLoadData, msg='Checking driver response type.')

            with driver_response as f:
                file_response_contents = f.read()

            self.assertEqual(reference_file_response_contents, file_response_contents,
                             msg='Checking file contents.')


class TestS3Driver(TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Set`s up mock S3 resource
        """
        # Init S3 mock
        cls.mock = mock_s3()
        cls.mock.start()

        cls.logger = logging.getLogger('s3_driver_test')
        cls.driver = S3Driver('tests3', 'tests/driver_configurations/aws_s3_config.yaml')

        # Add contents to mock S3 resource using boto3
        s3 = boto3.resource('s3')
        # Moto does not allow underscores in bucket names
        s3.create_bucket(Bucket='tests3', CreateBucketConfiguration={'LocationConstraint': 'us-west-2'})
        s3.Object('tests3', 'projects/').put(Body=b'')
        s3.Object('tests3', 'projects/test_file.txt').put(Body=b'S3 file projection testing data.',
                                                          Metadata={'madefor': 'testing', 'quality': 'good'})

    @classmethod
    def tearDownClass(cls):
        # Stopping mock resource
        cls.mock.stop()

    def test_bucket_contents(self):
        """
        Tests if driver gets proper bucket on init
        """
        exp_contents = ['projects/', 'projects/test_file.txt']

        self.assertEqual(len(exp_contents), len([o.key for o in self.driver.bucket.objects.all()]),
                         msg='Checking if none additional objects reported by driver')

        for obj in exp_contents:
            self.assertIn(obj, [o.key for o in self.driver.bucket.objects.all()], msg='Checking objects existence.')

    def test_get_uri_contents_as_dict(self):
        """
        Tests driver get_uri_contents_as_dict method
        """
        contents_meta = [{'content_encoding': None, 'size': 0, 'name': 'projects/',
                          'content_type': 'text/plain; charset=utf-8', 'metadata': {}, 'resource_uri': 'projects/'},
                         {'content_encoding': None, 'size': 32, 'name': 'projects/test_file.txt',
                          'content_type': 'text/plain; charset=utf-8',
                          'resource_uri': 'projects/test_file.txt',
                          'metadata': {'madefor': 'testing', 'quality': 'good'}}]
        for meta in contents_meta:
            self.assertDictEqual(meta, self.driver.get_uri_contents_as_dict(meta['name']),
                                 msg='Checking meta of object with URI: {0}'.format(meta['name']))

    def test_get_uri_contents_as_bytes(self):
        """
        Tests driver get_uri_contents_bytes method
        """
        test_data_contents = b'S3 file projection testing data.'
        for key, exp_content in {'projects/': b'', 'projects/test_file.txt': test_data_contents}.items():
            driver_response = self.driver.get_uri_contents_as_bytes(key)

            from drivers.aws_s3_driver import S3DataLoader
            self.assertIsInstance(driver_response, S3DataLoader, msg='Checking driver response type.')

            with driver_response as response:
                response = response.read()

            self.assertEqual(exp_content, response,
                             msg='Checking content of object with URI: {0}'.format(key))


class TestTorrentSuiteDriver(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_torrent_suite_driver')

        cls.mock_resource = MockResource('tests/torrent_suite_mock.json')

        cls.driver = TorrentSuiteDriver('mockiontorrent.com',
                                        'tests/driver_configurations/iontorrent_config.yaml')
        with open('tests/torrent_suite_mock.json') as m_f:
            cls.mock_structure = json.load(m_f)

        # Opening reference driver responses
        with open('tests/ts_driver_responses.json') as ts_dr_res:
            cls.reference_driver_responses = json.load(ts_dr_res)

    @classmethod
    def tearDownClass(cls):
        cls.mock_resource.deactivate()

    def test_get_uri_contents_as_dict(self):
        """
        Tests if torrent suite driver returns valid data on uri
        """
        for uri in self.reference_driver_responses:
            driver_response = self.driver.get_uri_contents_as_dict(uri)
            self.assertIsInstance(driver_response, dict, msg='Checking response object type.')

            self.assertDictEqual(self.reference_driver_responses[uri], driver_response,
                                 msg='Checking if current response is equal to example response')

    def test_get_uri_contents_as_bytes(self):
        """
        Tests if torrent suite driver gets contents of uri properly
        """
        uri_to_file_mapping = self.mock_structure['mock_responses']

        for uri in self.reference_driver_responses:
            driver_response = self.driver.get_uri_contents_as_bytes(uri)

            from drivers.iontorrent_driver import LoadTorrentSuiteData

            self.assertIsInstance(driver_response, LoadTorrentSuiteData,
                                  msg='Checking driver response type.')

            with driver_response as response:
                response = b''.join([el for el in response])

            for key in uri_to_file_mapping:
                # Getting corresponding data file for uri from mock structure
                if re.search(key, uri):
                    path_to_mock_contents_on_uri = os.path.join(self.mock_structure['mock_data_path'],
                                                                uri_to_file_mapping[key])

                    with open(path_to_mock_contents_on_uri, 'rb') as m_c:
                        mock_contents_on_uri = m_c.read()
                    self.assertEqual(mock_contents_on_uri, response, msg='Checking contents on uri.')


class TestGenbankDriver(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_torrent_suite_driver')

        cls.mock_resource = MockResource('tests/genbank_mock.json')

        cls.driver = GenbankDriver('http://eutils.ncbi.nlm.nih.gov',
                                   'tests/driver_configurations/genbank_config.yaml')

        with open('tests/genbank_driver_responses.json') as r_r:
            cls.reference_genbank_responses = json.load(r_r)

    @classmethod
    def tearDownClass(cls):
        cls.mock_resource.deactivate()

    def test_get_uri_contents_as_dict(self):
        """
        Tests if driver correctly returns dict on uri
        """
        for uri in self.reference_genbank_responses:
            driver_response = self.driver.get_uri_contents_as_dict(uri)

            self.assertIsInstance(driver_response, dict, msg='Checking type of driver response.')
            self.assertDictEqual(self.reference_genbank_responses[uri],
                                 driver_response,
                                 msg='Checking response content on uri:{}.'.format(uri))

    def test_get_uri_contents_as_bytes(self):
        """
        Tests if driver correctly returns uri contents
        """
        for uri in ['939732440.fasta', '939732440.json']:
            driver_response = self.driver.get_uri_contents_as_bytes(uri)

            if uri.endswith('.json'):
                reference_path = 'tests/mock_resources/genbank_mock_data/mock_json.json'
                from io import BytesIO
                self.assertIsInstance(driver_response, BytesIO, msg='Checking response type.')
                with driver_response as response:
                    response = b''.join([el for el in response])
                    response = json.loads(response.decode())
                with open(reference_path, 'r') as r_f:
                    reference_contents = json.load(r_f)
                self.assertEqual(reference_contents, response, msg='Checking contents on uri: {}.'.format(uri))

            elif uri.endswith('.fasta'):
                reference_path = 'tests/mock_resources/genbank_mock_data/mock_contents.txt'
                from drivers.genbank_driver import GenbankDataLoader
                self.assertIsInstance(driver_response, GenbankDataLoader, msg='Checking response type.')
                with open(reference_path, 'rb') as r_f:
                    reference_contents = r_f.read()
                logging.debug('Checking uri: %s', uri)
                with driver_response as response:
                    response = b''.join([el for el in response])
                self.assertEqual(reference_contents, response, msg='Checking contents on uri: {}.'.format(uri))


@skip("Driver creation is on hold due to sam-dump data loading difficulties.")
class TestSRADriver(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_torrent_suite_driver')

        cls.mock_resource = MockResource('tests/sra_mock.json')

        cls.driver = SRADriver('http://eutils.ncbi.nlm.nih.gov',
                               'tests/driver_configurations/sra_config.yaml')

        with open('tests/sra_driver_responses.json') as r_r:
            cls.reference_sra_responses = json.load(r_r)

    @classmethod
    def tearDownClass(cls):
        cls.mock_resource.deactivate()

    def test_get_uri_contents_as_dict(self):
        """
        Tests if driver correctly returns dict on uri
        """
        for uri in self.reference_sra_responses:
            driver_response = self.driver.get_uri_contents_as_dict(uri)

            self.assertIsInstance(driver_response, dict, msg='Checking type of driver response.')
            self.assertDictEqual(self.reference_sra_responses[uri],
                                 driver_response,
                                 msg='Checking response content.')
