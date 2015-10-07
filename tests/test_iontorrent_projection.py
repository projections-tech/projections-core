__author__ = 'abragin'

import os
import re
import logging
import logging.config
from unittest import TestCase, skip
from .mock import TorrentSuiteMock



import iontorrent

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('test_iontorrent_projection')

USER = 'user'
PASSWORD = 'password'


class TestIonTorrentProjection(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mock_resource = TorrentSuiteMock('mockiontorrent.com', 'tests/mock_resource')

    def setUp(self):
        self.mock_url = self.mock_resource.mock_url
        self.iontorrent = iontorrent.IonTorrentProjection(self.mock_url, USER, PASSWORD)
        self.iontorrent.create_projections()

    def test_bed_file_projections_creation(self):
        """
        Test BED file projection creation
        """
        bed_files_list = [p for p in self.iontorrent.projections if os.path.splitext(p)[1] == '.bed']
        expected_bed_file_path = os.path.join(self.expected_experiment_path,
                                              self.expected_run_name,
                                              'IAD39777_BED_4_for_TSVC.bed')

        self.assertTrue(expected_bed_file_path in bed_files_list)

    def test_full_projection(self):
        projection_paths_list = self.iontorrent.projections.keys()
        exp_dirs = ['test_experiment_1', 'test_experiment_2']
        run_name = 'test_run'
        # Checking number of created projections, we expect 372 projections for two experiments with 5 samples and variant_calling
        self.assertEqual(len(self.iontorrent.projections), 372, msg='Checking projections quantity.')

        for exp_dir in exp_dirs:
            # Checking BAM file projections creation on expected paths
            for i in range(1,6):
                bam_file_path = os.path.join('/', exp_dir, run_name, 'sample_{0}'.format(i), 'sample_{0}.bam'.format(i))
                self.assertTrue(bam_file_path in projection_paths_list, msg='Checking BAM projection existance: {0}'.format(bam_file_path))
                # Checking VCF files projection creation on expected path
                for variant_caller_dir in ['', '.49', '.50']:
                    for variant_file_name in ['TSVC_variants.vcf', 'all.merged.vcf', 'indel_assembly.vcf',
                                                              'indel_variants.vcf', 'small_variants.left.vcf',
                                                              'small_variants.vcf', 'small_variants_filtered.vcf',
                                                              'small_variants.sorted.vcf', 'SNP_variants.vcf']:
                        vcf_file_path = os.path.join('/', exp_dir,
                                                     run_name,
                                                     'sample_{}'.format(i),
                                                     'variantCaller_out{}'.format(variant_caller_dir),
                                                     variant_file_name)
                        self.assertTrue(vcf_file_path in projection_paths_list, msg='VCF path: {0}'.format(vcf_file_path))
        # Checking VCF files projections creation on expected paths


