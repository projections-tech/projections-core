__author__ = 'abragin'

import os
import logging
import logging.config
from unittest import TestCase, skip
from tests.torrent_suite_mock import TorrentSuiteMock
from projections import PrototypeDeserializer, Projection

import iontorrent


MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_iontorrent_projection')

USER = 'user'
PASSWORD = 'password'


class TestTorrentSuiteProjector(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mock_resource = TorrentSuiteMock('mockiontorrent.com', 'tests/mock_resource/torrent_suite_mock_data')

    def test_full_projection(self):
        """
        Tests projections correctness of full projection from root creation.
        """
        projection_configuration = PrototypeDeserializer('tests/test_full_torrent_suite_config.yaml')
        root_projection = Projection('/', projection_configuration.root_projection_uri)
        driver = iontorrent.TorrentSuiteDriver(projection_configuration.resource_uri, USER, PASSWORD)
        iontorrent_projector = iontorrent.TorrentSuiteProjector(driver, root_projection,
                                                           projection_configuration.prototype_tree)

        projection_paths_list = iontorrent_projector.projections.keys()

        # Checking number of created projections,
        self.assertEqual(len(projection_paths_list), 131,
                         msg='Checking total number of projections,'
                             ' expecting 131, got: {}.'.format(len(projection_paths_list)))

        # Representation of mock internal structure
        experiment_contents = {
            '/test_experiment_1':
                {
                    'test_run_1': ['sample_1', 'sample_2'],
                    'test_run_2': ['sample_1', 'sample_2']
                },
            '/test_experiment_2':
                {
                    'test_run_3': ['sample_3', 'sample_4'],
                    'test_run_4': ['sample_3', 'sample_4']
                }
        }

        for exp_dir in self.mock_resource.get_experiments():
            exp_dir = '/' + exp_dir
            # Checking metadata projection creation for experiments
            for meta_name in ['metadata.json', 'plannedexperiment.json']:
                meta_data_path = os.path.join(exp_dir, meta_name)
                self.assertTrue(meta_data_path in projection_paths_list,
                                msg='Checking metadata projections for experiment {}.'.format(exp_dir))

            for run_name, sample_ids in experiment_contents[exp_dir].items():
                for sample_id in sample_ids:
                    # Checking BAM file projections creation on expected paths
                    bam_file_path = os.path.join(exp_dir, run_name, sample_id, '{0}.bam'.format(sample_id))
                    self.assertIn(bam_file_path, projection_paths_list,
                                  msg='Checking BAM projection existence: {}'.format(bam_file_path))

                    # Test sample metadata projection creation
                    sample_meta_path = os.path.join(exp_dir, run_name, sample_id, 'metadata.json')
                    self.assertIn(sample_meta_path, projection_paths_list,
                                  msg='Checking metadata creation:{}'.format(sample_meta_path))

                    variant_caller_dir = '.587'
                    vc_dir_path = os.path.join(exp_dir, run_name, sample_id,
                                               'variantCaller_out{}'.format(variant_caller_dir))

                    # Checking BED file projections creation
                    bed_file_path = os.path.join(vc_dir_path, 'IAD66589_181_SNP-HID-p2-L_Target_regions.bed')
                    self.assertIn(bed_file_path, projection_paths_list,
                                  msg='Checking BED file projection existence: {}.'.format(bed_file_path))

                    # Checking variant caller settings projection creation
                    vcf_settings_path = os.path.join(vc_dir_path, 'variant_caller_settings.json')
                    self.assertIn(vcf_settings_path, projection_paths_list,
                                  msg='Checking VC settings projection creation: {}'.format(vcf_settings_path))

                    # Checking VCF files projection creation on expected path
                    for variant_file_name in ['TSVC_variants.vcf', 'all.merged.vcf', 'indel_assembly.vcf',
                                              'indel_variants.vcf', 'small_variants.left.vcf',
                                              'small_variants.vcf', 'small_variants_filtered.vcf',
                                              'small_variants.sorted.vcf', 'SNP_variants.vcf']:
                        vcf_file_path = os.path.join(exp_dir,
                                                     run_name,
                                                     sample_id,
                                                     'variantCaller_out{}'.format(variant_caller_dir),
                                                     variant_file_name)
                        self.assertIn(vcf_file_path, projection_paths_list, msg='VCF path: {}'.format(vcf_file_path))

    def test_non_root_projection(self):
        """
        Test creation of projection with user specified root
        """

        # Loading configuration where '/rundb/api/v1/results/1/' is root projection
        projection_configuration = PrototypeDeserializer('tests/test_custom_root_torrent_suite_config.yaml')
        root_projection = Projection('/', projection_configuration.root_projection_uri)
        driver = iontorrent.TorrentSuiteDriver(projection_configuration.resource_uri, USER, PASSWORD)
        iontorrent_projector = iontorrent.TorrentSuiteProjector(driver, root_projection,
                                                                projection_configuration.prototype_tree)

        projection_paths_list = iontorrent_projector.projections.keys()

        # Checking number of created projections,
        self.assertEqual(len(projection_paths_list), 32,
                         msg='Checking total number of projections,'
                             ' expecting 32, got: {}.'.format(len(projection_paths_list)))

        run_name = '/test_run_1'
        sample_ids = ['sample_1', 'sample_2']

        for sample_id in sample_ids:
            # Checking BAM file projections creation on expected paths
            bam_file_path = os.path.join(run_name, sample_id, '{0}.bam'.format(sample_id))
            self.assertIn(bam_file_path, projection_paths_list,
                          msg='Checking BAM projection existence: {}'.format(bam_file_path))

            # Test sample metadata projection creation
            sample_meta_path = os.path.join(run_name, sample_id, 'metadata.json')
            self.assertIn(sample_meta_path, projection_paths_list,
                          msg='Checking metadata creation:{}'.format(sample_meta_path))

            variant_caller_dir = '.587'
            vc_dir_path = os.path.join(run_name, sample_id,
                                       'variantCaller_out{}'.format(variant_caller_dir))

            # Checking BED file projections creation
            bed_file_path = os.path.join(vc_dir_path, 'IAD66589_181_SNP-HID-p2-L_Target_regions.bed')
            self.assertIn(bed_file_path, projection_paths_list,
                          msg='Checking BED file projection existence: {}.'.format(bed_file_path))

            # Checking variant caller settings projection creation
            vcf_settings_path = os.path.join(vc_dir_path, 'variant_caller_settings.json')
            self.assertIn(vcf_settings_path, projection_paths_list,
                          msg='Checking VC settings projection creation: {}'.format(vcf_settings_path))

            # Checking VCF files projection creation on expected path
            for variant_file_name in ['TSVC_variants.vcf', 'all.merged.vcf', 'indel_assembly.vcf',
                                      'indel_variants.vcf', 'small_variants.left.vcf',
                                      'small_variants.vcf', 'small_variants_filtered.vcf',
                                      'small_variants.sorted.vcf', 'SNP_variants.vcf']:
                vcf_file_path = os.path.join(run_name,
                                             sample_id,
                                             'variantCaller_out{}'.format(variant_caller_dir),
                                             variant_file_name)
                self.assertIn(vcf_file_path, projection_paths_list, msg='VCF path: {}'.format(vcf_file_path))

    def test_projection_contraction(self):
        projection_configuration = PrototypeDeserializer('tests/test_torrent_suite_projection_filtering_config.yaml')
        root_projection = Projection('/', projection_configuration.root_projection_uri)
        driver = iontorrent.TorrentSuiteDriver(projection_configuration.resource_uri, USER, PASSWORD)
        iontorrent_projector = iontorrent.TorrentSuiteProjector(driver, root_projection,
                                                           projection_configuration.prototype_tree)

        projection_paths_list = iontorrent_projector.projections.keys()


        # Checking number of created projections,
        self.assertEqual(len(projection_paths_list), 19,
                         msg='Checking total number of projections,'
                             ' expecting 19, got: {}.'.format(len(projection_paths_list)))

        # Representation of mock internal structure
        experiment_contents = {
                    'test_run_1': ['sample_1', 'sample_2']
                }

        exp_dir = '/test_experiment_1'

        self.assertIn(exp_dir, projection_paths_list, msg='Checking experiment 1 dir projection creation')

        # Checking experiment filtration by name
        self.assertNotIn('/test_experiment_2', projection_paths_list, msg='Checking experiment 2 dir projection is not created')

        # Checking run filtration by name
        self.assertNotIn('/test_experiment_1/test_run_2', projection_paths_list, msg='Checking test run 2 dir projection is not created')

        # Checking metadata projection creation for experiments
        for meta_name in ['metadata.json', 'plannedexperiment.json']:
            meta_data_path = os.path.join(exp_dir, meta_name)
            self.assertTrue(meta_data_path in projection_paths_list,
                            msg='Checking metadata projections for experiment {}.'.format(exp_dir))

        for run_name, sample_ids in experiment_contents.items():
            for sample_id in sample_ids:
                # Checking BAM file projections creation on expected paths
                bam_file_path = os.path.join(exp_dir, run_name, sample_id, '{0}.bam'.format(sample_id))
                self.assertIn(bam_file_path, projection_paths_list,
                              msg='Checking BAM projection existence: {}'.format(bam_file_path))

                # Test sample metadata projection creation
                sample_meta_path = os.path.join(exp_dir, run_name, sample_id, 'metadata.json')
                self.assertIn(sample_meta_path, projection_paths_list,
                              msg='Checking metadata creation:{}'.format(sample_meta_path))

                variant_caller_dir = '.587'
                vc_dir_path = os.path.join(exp_dir, run_name, sample_id,
                                           'variantCaller_out{}'.format(variant_caller_dir))

                # Checking BED file projections creation
                bed_file_path = os.path.join(vc_dir_path, 'IAD66589_181_SNP-HID-p2-L_Target_regions.bed')
                self.assertIn(bed_file_path, projection_paths_list,
                              msg='Checking BED file projection existence: {}.'.format(bed_file_path))

                # Checking variant caller settings projection creation
                vcf_settings_path = os.path.join(vc_dir_path, 'variant_caller_settings.json')
                self.assertIn(vcf_settings_path, projection_paths_list,
                              msg='Checking VC settings projection creation: {}'.format(vcf_settings_path))

                # Checking if selected TSVC VCF file was created
                tsvc_vcf_file_path = os.path.join(exp_dir, run_name, sample_id,
                                             'variantCaller_out{}'.format(variant_caller_dir),
                                             'TSVC_variants.vcf')
                self.assertIn(tsvc_vcf_file_path, projection_paths_list, msg='VCF path: {}'.format(tsvc_vcf_file_path))

                # Checking VCF files filtration
                for variant_file_name in ['all.merged.vcf', 'indel_assembly.vcf',
                                          'indel_variants.vcf', 'small_variants.left.vcf',
                                          'small_variants.vcf', 'small_variants_filtered.vcf',
                                          'small_variants.sorted.vcf', 'SNP_variants.vcf']:
                    vcf_file_path = os.path.join(exp_dir,
                                                 run_name,
                                                 sample_id,
                                                 'variantCaller_out{}'.format(variant_caller_dir),
                                                 variant_file_name)
                    self.assertNotIn(vcf_file_path, projection_paths_list, msg='VCF path: {}'.format(vcf_file_path))