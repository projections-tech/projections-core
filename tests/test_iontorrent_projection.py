__author__ = 'abragin'

import os
import logging
import logging.config
from unittest import TestCase, skip
from tests.torrent_suite_mock import TorrentSuiteMock
from projections import PrototypeDeserializer

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

    def setUp(self):
        projection_configuration = PrototypeDeserializer('torrent_suite_config.yaml')
        driver = iontorrent.TorrentSuiteDriver(projection_configuration.resource_uri, USER, PASSWORD)
        self.iontorrent = iontorrent.TorrentSuiteProjector(driver, projection_configuration.prototype_tree)

    def test_full_projection(self):
        """
        Tests projections correctness of full projection creation.
        """
        projection_paths_list = self.iontorrent.projections.keys()
        exp_dirs = ['/test_experiment_1', '/test_experiment_2']
        run_name = 'test_run'

        # Checking number of created projections,
        # we expect 39 projections for two experiments with 1 sample and variant_calling
        self.assertEqual(len(projection_paths_list), 39,
                         msg='Checking total number of projections,'
                             ' expecting 39, got: {}.'.format(len(projection_paths_list)))

        for exp_dir in exp_dirs:
            # Checking metadata projection creation for experiments
            for meta_name in ['metadata.json', 'plannedexperiment.json']:
                meta_data_path = os.path.join(exp_dir, meta_name)
                self.assertTrue(meta_data_path in projection_paths_list,
                                msg='Checking metadata projections for experiment {}.'.format(exp_dir))

            # Checking BAM file projections creation on expected paths
            bam_file_path = os.path.join(exp_dir, run_name, 'sample_1', 'sample_1.bam')
            self.assertIn(bam_file_path, projection_paths_list,
                          msg='Checking BAM projection existence: {}'.format(bam_file_path))

            # Test sample metadata projection creation
            sample_meta_path = os.path.join(exp_dir, run_name, 'sample_1', 'metadata.json')
            self.assertIn(sample_meta_path, projection_paths_list,
                          msg='Checking metadata creation:{}'.format(sample_meta_path))

            variant_caller_dir = '.587'
            vc_dir_path = os.path.join(exp_dir, run_name, 'sample_1',
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
                                             'sample_1',
                                             'variantCaller_out{}'.format(variant_caller_dir),
                                             variant_file_name)
                self.assertIn(vcf_file_path, projection_paths_list, msg='VCF path: {}'.format(vcf_file_path))