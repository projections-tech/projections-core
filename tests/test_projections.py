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

import getpass
import json
import logging
import logging.config
import re
from io import BytesIO
from unittest import TestCase

import psycopg2

from db_projector import DBProjector
from drivers.fs_driver import FSDriver
from projections import ProjectionPrototype, ProjectionDriver, PrototypeDeserializer

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')


# TODO implement this class as mock from unittest library
class TestDriver(ProjectionDriver):
    """
    This class does not perform actual testing. It's providing test data and may be replaced with mock object.
    """

    def __init__(self):
        self.logger = logging.getLogger('test_driver')

    def get_uri_contents_as_dict(self, uri):
        self.logger.debug('Requesting content for uri: %s', uri)

        if not uri:
            self.logger.debug('Returning empty array')
            return []

        if uri == 'experiments':
            with open('tests/json/experiments.json') as f:
                content = json.load(f)
            self.logger.debug('Returning content for uri: %s, content: %s', uri, content)
            return content

        match = re.match('experiments/(\d+)', uri)
        if match:
            id = match.groups()[0]
            self.logger.debug('Requesting experiment data with id: %s', id)
            with open('tests/json/experiment_{}.json'.format(id)) as f:
                content = json.load(f)
            self.logger.debug('Returning content for uri: %s, content: %s', uri, content)
            return content

        match = re.match('results/(\d+)', uri)
        if match:
            id = match.groups()[0]
            self.logger.debug('Requesting result data with id: %s', id)
            return {'id': id, 'content': 'Result of some experiment',
                    'filesystempath': '/tmp/result_{}'.format(id), 'data': 'data/{}.bam'.format(id)}

        match = re.match('data/(\d+).bam', uri)
        if match:
            id = match.groups()[0]
            self.logger.debug('Requesting result data with id: %s', id)
            return {'meta': "This is BAM file"}

        assert False is True, 'Test driver can\'t handle resource request, aborting!'

    def get_uri_contents_as_bytes(self, uri):
        return TestDriverLoadData(uri)


class TestDriverLoadData:
    """
    Data loading context manager for test driver
    """

    def __init__(self, uri):
        self.uri = uri
        self.io = None

    def __enter__(self):
        self.io = BytesIO(b'Mock resource contents.')
        return self.io

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.io.close()

class TestPrototypeDeserializer(TestCase):
    def setUp(self):
        self.deserializer = PrototypeDeserializer('tests/test_projection_config.yaml')

    def test_prototype_deserialization(self):
        """
        Tests test_projection_config deserialization in to ProjectionPrototype tree
        """
        root_prototype = self.deserializer.prototype_tree
        # Test if nodes are ProjectionPrototype instances
        test_nodes_list = [n for n in root_prototype.get_tree_nodes()]
        for n in test_nodes_list:
            self.assertIsInstance(n, ProjectionPrototype,
                                  msg='Checking if object: {0} is instance of ProjectionPrototype'.format(
                                      root_prototype))
        # Test correctness of "name" fields of nodes
        expected_names = ['root_dir', 'results_dir', 'test_bam.bam', 'test_vcf.vcf']
        test_names = [n.name for n in root_prototype.get_tree_nodes()]
        for element in expected_names:
            self.assertIn(element, test_names,
                          msg='Checking existance of projection with name {} in a tree.'.format(element))
        # Test correctness of "uri" fields of nodes
        expected_uri = ['[object["uri"] for object in environment]', "environment['results']",
                        "[environment['data_vcf']]", "[environment['data_bam']]"]
        test_pre_order_uri = [n.uri for n in root_prototype.get_tree_nodes()]
        for element in expected_uri:
            self.assertIn(element, test_pre_order_uri,
                          msg='Checking existance of projection with uri {} in a tree.'.format(element))
        # Test correctness of "type" fields of nodes
        expected_types = ['directory', 'directory', 'file', 'file']
        test_pre_order_uri = [n.type for n in root_prototype.get_tree_nodes()]
        self.assertListEqual(expected_types, test_pre_order_uri, msg='Checking if prototypes types are correct.')


class TestProjector(TestCase):
    """
    Tests DBProjector projection tree building
    """

    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_projector')

        # Initializing database connection which will be used during tests
        cls.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))
        # Creating cursor, which will be used to interact with database
        cls.cursor = cls.db_connection.cursor()

        cls.cursor.execute(" DELETE FROM projections.projections WHERE projection_name='test_projection'; ")
        cls.db_connection.commit()

        cls.projection_driver = TestDriver()

    @classmethod
    def tearDownClass(cls):
        cls.cursor.execute(" DELETE FROM projections.projections WHERE projection_name='test_projection'; ")
        cls.db_connection.commit()

        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        self.cursor.execute("""
        INSERT INTO projections.projections (projection_name, mount_point, driver, prototype)
        VALUES ('test_projection', 'None', 'test_driver', 'test_prototype')
        RETURNING projection_id
        """)

        self.projection_id = self.cursor.fetchone()
        self.db_connection.commit()

        projection_settings = PrototypeDeserializer('tests/projections_configs/test_projection.yaml')

        self.projector = DBProjector(self.projection_id,
                                     self.projection_driver,
                                     projection_settings.prototype_tree,
                                     projection_settings.root_projection_uri)

    def tearDown(self):
        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM projections.projections WHERE projection_name='test_projection'; ")
        self.db_connection.commit()

    def _list_projections(self):
        """
        This function is used to return list of projections paths from tree table
        :returns: list of strings
        """
        self.cursor.execute("""
        SELECT join_path(node_path, node_name)
        FROM projections.projection_nodes
        WHERE tree_id=%s
        """, (self.projection_id))

        return [row[0] for row in self.cursor]

    def test_db_build_tree(self):
        """
        Testing projection tree creation with projection prototypes.
        """

        dir_paths = ['/', '/experiment_0', '/experiment_1', '/experiment_2',
                     '/experiment_1/result_1', '/experiment_1/result_2',
                     '/experiment_2/result_3', '/experiment_2/result_4', '/experiment_2/result_5']

        created_projections = self._list_projections()

        self.logger.debug('Created projections: %s', created_projections)
        for dir_path in dir_paths:
            self.logger.debug('Checking projection on path: %s', dir_path)

            self.assertIn(dir_path, created_projections, 'Check that projection exists')

            projection_stats = self.projector.get_attributes(dir_path)
            self.logger.debug('Proj stats: %s', projection_stats)
            self.assertEqual(projection_stats['st_mode'], 16895, 'Check that this is a directory projection')

        file_paths = ['/experiment_1/result_1/1.bam', '/experiment_1/result_2/2.bam',
                      '/experiment_2/result_3/3.bam', '/experiment_2/result_4/4.bam', '/experiment_2/result_5/5.bam']

        for file_path in file_paths:
            self.logger.debug('Checking file projection on path: %s', file_path)
            self.assertIn(file_path, created_projections, 'Check that projection exists')

            projection_stats = self.projector.get_attributes(file_path)

            self.assertTrue(projection_stats['st_mode'] == 33279, 'Check that this is a file projection')

    def test_remove_projection(self):
        """
        Tests projection deletion correctness
        """
        # Removing entire directory
        self.projector.remove_projection('/experiment_1')

        current_projections = self._list_projections()

        removed_projections = ['/experiment_1', '/experiment_1/result_1',
                               '/experiment_1/result_1/1.bam', '/experiment_1/result_2',
                               '/experiment_1/result_2/2.bam']

        # Checking if projections where removed properly
        for removed_projection in removed_projections:
            self.assertNotIn(removed_projection, current_projections, 'Checking directory deletion')

        # Testing leaf projections removal
        projections_to_remove = ['/experiment_2/result_4/4.bam', '/experiment_2/result_5/5.bam']
        # Removing leaf projections
        for path in projections_to_remove:
            self.projector.remove_projection(path)

        current_projections = self._list_projections()

        for path in projections_to_remove:
            self.assertNotIn(path, current_projections, 'Checking file deletion')

    def test_move_projection(self):
        """
        Tests projections move correctness
        """
        # Moving dir experiment_1 into dir experiment_2
        self.projector.move_projection('/experiment_1', '/experiment_2')

        current_projections = self._list_projections()

        moved_dir_paths = ['/experiment_2/experiment_1', '/experiment_2/experiment_1/result_2',
                           '/experiment_2/experiment_1/result_1',
                           '/experiment_2/experiment_1/result_2/2.bam',
                           '/experiment_2/experiment_1/result_1/1.bam']
        for path in moved_dir_paths:
            self.assertIn(path, current_projections, 'Checing node move')

    def test_is_managing_path(self):
        """
        Test if projection manager reports projection managment status properly.
        """
        projection_paths = ['/', '/experiment_0', '/experiment_1', '/experiment_2',
                            '/experiment_1/result_1', '/experiment_1/result_2',
                            '/experiment_2/result_3', '/experiment_2/result_4',
                            '/experiment_2/result_5', '/experiment_1/result_1/1.bam',
                            '/experiment_1/result_2/2.bam', '/experiment_2/result_3/3.bam',
                            '/experiment_2/result_4/4.bam', '/experiment_2/result_5/5.bam']
        for path in projection_paths:
            self.assertTrue(self.projector.is_managing_path(path), msg='Testing if projector manages path.')

    def test_get_attributes(self):
        """
        Tests if projector reports correct projections attributes
        """
        created_projections = self._list_projections()
        for path in created_projections:
            self.cursor.execute("""
            SELECT st_atime, st_mtime, st_ctime, st_size, st_mode, st_nlink, st_ino
            FROM projections.get_projection_node_attributes(%s, %s)
            """, (self.projection_id, path))

            for row in self.cursor:
                st_atime, st_mtime, st_ctime, st_size, st_mode, st_nlink, st_ino = row
                # Checking attributes type correctness
                self.assertIsInstance(st_atime, int, msg='Checking projection st_atime attribute type.')
                self.assertIsInstance(st_mtime, int, msg='Checking projection st_mtime attribute type.')
                self.assertIsInstance(st_ctime, int, msg='Checking projection st_ctime attribute type.')
                self.assertIsInstance(st_size, int, msg='Checking projection st_size attribute type.')
                self.assertIsInstance(st_nlink, int, msg='Checking projection st_nlink attribute type.')
                self.assertIsInstance(st_ino, int, msg='Checking projection st_ino attribute type.')
                self.assertIsInstance(st_mode, str, msg='Checking projection st_mode attribute type.')
                # Checking if attributes returned correctly
                self.assertEqual(st_size, 1, msg='Checkin projection st_size attribute value.')
                self.assertEqual(st_nlink, 0, msg='Checkin projection st_nlink attribute value.')
                self.assertEqual(st_ino, 1, msg='Checkin projection st_ino attribute value.')

    def test_get_projections_on_path(self):
        """
        Tests if projector correctly returns projections on path
        """
        self.assertEqual(set(['experiment_0', 'experiment_1', 'experiment_2']),
                         set(self.projector.get_projections_on_path('/')),
                         msg='Checking get_projections_on_path on dir path.')

        self.assertListEqual([], self.projector.get_projections_on_path('/experiment_1/result_2/2.bam'),
                             msg='Checking get_projections_on_path on file path.')

    def test_open_resource(self):
        """
        Tests if projector opens resources properly
        """
        paths_list = self._list_projections()
        for path in paths_list:
            header, context_manager = self.projector.open_resource(path)

            self.assertIsInstance(context_manager, TestDriverLoadData,
                                  msg='Checks if projector returns driver context manager.')

            with context_manager as manager:
                path_data = manager.read()

            self.assertEqual(b'Mock resource contents.', path_data, msg='Check resource contents.')

            self.assertIsInstance(header, int, msg='Check header type.')
            self.assertEqual(3, header, msg='Check header value.')


class TestProjectionVariants(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_projector')

        # Initializing database connection which will be used during tests
        cls.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))
        # Creating cursor, which will be used to interact with database
        cls.cursor = cls.db_connection.cursor()

        cls.cursor.execute(" DELETE FROM projections.projections WHERE projection_name='test_projection'; ")
        cls.db_connection.commit()

        cls.projection_driver = TestDriver()

    @classmethod
    def tearDownClass(cls):
        cls.cursor.execute(" DELETE FROM projections.projections WHERE projection_name='test_projection'; ")
        cls.db_connection.commit()

        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        self.cursor.execute("""
        INSERT INTO projections.projections (projection_name, mount_point, driver, prototype)
        VALUES ('test_projection', 'None', 'test_driver', 'test_prototype')
        RETURNING projection_id
        """)
        self.projection_id = self.cursor.fetchone()
        self.db_connection.commit()

    def tearDown(self):
        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM projections.projections WHERE projection_name='test_projection'; ")
        self.db_connection.commit()

    def _list_projections(self):
        """
        This function is used to return list of projections paths from tree table
        :returns: list of strings
        """
        self.cursor.execute("""
        SELECT join_path(node_path, node_name)
        FROM projections.projection_nodes
        WHERE tree_id=%s
        """, (self.projection_id,))

        return [row[0] for row in self.cursor]

    def test_transparent_projection_creation(self):
        """
        Tests use of transparent projections
        """

        projection_settings = PrototypeDeserializer('tests/projections_configs/test_transaprent_projection_config.yaml')

        projector = DBProjector(self.projection_id,
                                self.projection_driver,
                                projection_settings.prototype_tree,
                                projection_settings.root_projection_uri)

        expected_projections = ['/', '/result_1_1.bam', '/result_2_2.bam',
                                '/result_3_3.bam', '/result_4_4.bam', '/result_5_5.bam']

        created_projections = self._list_projections()

        self.assertEqual(set(expected_projections), set(created_projections),
                         msg='Checking transparent projections creation.')

    def test_projection_filtration(self):
        """
        Tests projection filtration on prototype microcode level
        """

        projection_settings = PrototypeDeserializer('tests/projections_configs/test_projection_filtration_config.yaml')

        projector = DBProjector(self.projection_id,
                                self.projection_driver,
                                projection_settings.prototype_tree,
                                projection_settings.root_projection_uri)

        expected_projections = ['/', '/experiment_1', '/experiment_1/result_1', '/experiment_1/result_1/1.bam',
                                '/experiment_1/result_2', '/experiment_1/result_2/2.bam']

        created_projections = self._list_projections()

        self.assertEqual(set(expected_projections), set(created_projections),
                         msg='Checking projection filtration.')

    def test_non_root_resource_projection(self):
        """
        Tests non root resource projection creation
        """
        projection_settings = PrototypeDeserializer('tests/projections_configs/test_non_root_projection.yaml')

        projector = DBProjector(self.projection_id,
                                self.projection_driver,
                                projection_settings.prototype_tree,
                                projection_settings.root_projection_uri)

        expected_projections = ['/', '/result_1', '/result_1/1.bam', '/result_2', '/result_2/2.bam']

        created_projections = self._list_projections()

        self.assertEqual(set(expected_projections), set(created_projections),
                         msg='Checking non root resource projection creation.')


class TestMetadataOperations(TestCase):
    """
    Tests DBProjector projection tree building
    """

    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_projector')

        # Initializing database connection which will be used during tests
        cls.db_connection = psycopg2.connect(
            "dbname=projections_database user={user_name}".format(user_name=getpass.getuser()))
        # Creating cursor, which will be used to interact with database
        cls.cursor = cls.db_connection.cursor()

        cls.cursor.execute(" DELETE FROM projections.projections WHERE projection_name='test_projection'; ")
        cls.db_connection.commit()

        cls.projection_driver = TestDriver()

    @classmethod
    def tearDownClass(cls):
        cls.cursor.execute(" DELETE FROM projections.projections WHERE projection_name='test_projection'; ")
        cls.db_connection.commit()

        # Closing cursor and connection
        cls.cursor.close()
        cls.db_connection.close()

    def setUp(self):
        self.cursor.execute("""
        INSERT INTO projections.projections (projection_name, mount_point, driver, prototype)
        VALUES ('test_projection', 'None', 'test_driver', 'test_prototype')
        RETURNING projection_id
        """)
        self.projection_id = self.cursor.fetchone()
        self.db_connection.commit()

        projection_settings = PrototypeDeserializer('tests/projections_configs/test_metadata_operations.yaml')

        projection_driver = FSDriver(projection_settings.root_projection_uri, projection_settings.driver_config_path)

        self.projector = DBProjector(self.projection_id,
                                     projection_driver,
                                     projection_settings.prototype_tree,
                                     projection_settings.root_projection_uri)

    def tearDown(self):
        # Clean up previous test entries in db
        self.cursor.execute(" DELETE FROM projections.projections WHERE projection_name='test_projection'; ")
        self.db_connection.commit()

    def _list_projections(self):
        """
        This function is used to return list of projections paths from tree table
        :returns: list of strings
        """
        self.cursor.execute("""
        SELECT join_path(node_path, node_name)
        FROM projections.projection_nodes
        WHERE tree_id=%s
        """, (self.projection_id,))

        return [row[0] for row in self.cursor]

    def test_metadata_binding(self):
        """
        Tests metadata-data binding process.
        """
        # Loading data-metadata name pairs
        self.cursor.execute("""
        WITH parents AS (
        SELECT projections.projection_links.tail_node_id as meta_id,
                projections.projection_links.head_node_id as parent_id,
                projections.projection_nodes.node_name as parent_name
        FROM projections.projection_nodes, projections.projection_links
        WHERE projections.projection_nodes.node_id = projections.projection_links.head_node_id
        )
        SELECT parents.parent_name, projections.projection_nodes.node_name
        FROM projections.projection_nodes, parents
        WHERE projections.projection_nodes.node_id = parents.meta_id
        """)

        parent_meta_pairs = [row for row in self.cursor]

        for row in parent_meta_pairs:
            parent_name, metadata_name = row
            self.logger.debug('Parent_name: %s Meta_name: %s', parent_name, metadata_name)
            if parent_name.endswith('.bam'):
                parent_name = parent_name.replace('.bam', '')
                metadata_name = metadata_name.replace('_metadata.json', '')
                self.assertEqual(parent_name, metadata_name, msg='Checking correctness of BAM metadata binding.')
            elif parent_name.endswith('.fasta'):
                self.assertEqual('fasta_file_.fasta', re.sub('\d+', '', parent_name),
                                 msg='Checking if parent node is fasta file.')
                self.assertEqual('vcf_file.vcf', metadata_name,
                                 msg='Checking if meta node is vcf file.')
