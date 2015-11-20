import os
import sys
import time
import subprocess
import logging
import logging.config
from unittest import TestCase, skip
from projections import PrototypeDeserializer, Projector
from fs_projection import FSDriver

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'
CONFIG_PATH = 'tests/test_fs_proj_config.yaml'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('fs_proj_test')


class TestFSProjection(TestCase):

    def setUp(self):
        driver = FSDriver()
        projection_configuration = PrototypeDeserializer(CONFIG_PATH)
        self.fs_projector = Projector(driver, projection_configuration.root_projection_uri,
                                      projection_configuration.prototype_tree)

    def test_create_projections(self):
        """
        Tests if filesystem projector creates projections
        """

        created_projections = [n.get_path() for n in self.fs_projector.projection_tree.get_tree_nodes()]

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(37, len(created_projections),
                         msg='FS projector created {0} projections.'.format(len(created_projections)))

        # Check root dir projection creation
        self.assertIn('/test_dir',
                      created_projections,
                      msg='Checking creation of root projection')

        # Checking if FASTA files not present in config are not projected
        for i in range(1,5):
            fasta_proj_path = '/test_dir/fasta_file_{0}.fasta'.format(i)
            self.assertNotIn(fasta_proj_path, created_projections,
                          msg='Checking if {0} not in projections'.format(fasta_proj_path))

        # Check BAM files projections creation
        for i in range(1,6):
            bam_proj_path = '/test_dir/bam_file_{0}.bam'.format(i)
            self.assertIn(bam_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(bam_proj_path))

        for i in range(1,6):
            # Check sample dir projection creation
            sample_dir_path = '/test_dir/sample_{0}'.format(i)
            self.assertIn(sample_dir_path, created_projections,
                          msg='Checking creation of {0} projection'.format(sample_dir_path))

            # Check sample VCF and BED projections creation
            vcf_proj_path = os.path.join(sample_dir_path, 'vcf_file.vcf')
            bed_proj_path = os.path.join(sample_dir_path, 'bed_file.bed')
            self.assertIn(vcf_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(vcf_proj_path))
            self.assertIn(bed_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(bed_proj_path))

            # Check inner rerun projection creation
            rerun_proj_path = os.path.join(sample_dir_path, 'rerun')
            self.assertIn(rerun_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(rerun_proj_path))

            # Check rerun contents BED and VCF files
            vcf_proj_path = os.path.join(rerun_proj_path, 'vcf_file.vcf')
            bed_proj_path = os.path.join(rerun_proj_path, 'bed_file.bed')
            self.assertIn(vcf_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(vcf_proj_path))
            self.assertIn(bed_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(bed_proj_path))

    def test_projections_contents(self):
        """
        Tests if created projections have proper content
        """
        # Unmount any previous tests
        subprocess.Popen(['fusermount', '-u', MOUNT_POINT])

        fs_proj = subprocess.Popen([sys.executable,
                                    'fs_projection.py',
                                    '-m', MOUNT_POINT,
                                    '-d', DATA_FOLDER,
                                    '-c', CONFIG_PATH],
                                    stdout=subprocess.DEVNULL)
        # Time to initialize Projector properly
        time.sleep(0.2)

        # Expected mock files contents
        expected_bam = 'Mock bam here!\n'
        expected_vcf = 'Mock vcf here!\n'
        expected_bed = 'Mock bed here!\n'
        for i in range(1,6):
            # Paths to projections
            vcf_proj_path = 'tests/mnt/test_dir/sample_{}/vcf_file.vcf'.format(i)
            bam_proj_path = 'tests/mnt/test_dir/bam_file_{0}.bam'.format(i)
            bed_proj_path = 'tests/mnt/test_dir/sample_{}/bed_file.bed'.format(i)
            rerun_vcf_proj_path = 'tests/mnt/test_dir/sample_{}/rerun/vcf_file.vcf'.format(i)
            rerun_bed_proj_path = 'tests/mnt/test_dir/sample_{}/rerun/bed_file.bed'.format(i)

            # Reading projections contents on paths
            with open(bam_proj_path) as bam_f:
                test_bam_file = bam_f.read()

            with open(vcf_proj_path) as vcf_f:
                test_vcf_file = vcf_f.read()

            with open(bed_proj_path) as bed_f:
                test_bed_file = bed_f.read()

            with open(rerun_vcf_proj_path) as r_vcf_f:
                test_r_vcf_file = r_vcf_f.read()

            with open(rerun_bed_proj_path) as r_bed_f:
                test_r_bed_file = r_bed_f.read()

            # Checking projected files contents
            self.assertEqual(expected_bam, test_bam_file,
                             msg='Checking contents of projection {0}'.format(bam_proj_path))
            self.assertEqual(expected_vcf, test_vcf_file,
                             msg='Checking contents of projection {0}'.format(vcf_proj_path))
            self.assertEqual(expected_bed, test_bed_file,
                             msg='Checking contents of projection {0}'.format(bed_proj_path))
            self.assertEqual(expected_vcf, test_r_vcf_file,
                             msg='Checking contents of projection {0}'.format(rerun_vcf_proj_path))
            self.assertEqual(expected_bed, test_r_bed_file,
                             msg='Checking contents of projection {0}'.format(rerun_bed_proj_path))

        # Terminating Projector process
        fs_proj.terminate()
        # Unmounting mnt dir
        subprocess.Popen(['fusermount', '-u', MOUNT_POINT])

    def test_projection_with_file_metadata(self):
        """
        Tests projection filtration according to data in metadata files placed with projected file
        """
        driver = FSDriver()
        projection_configuration = PrototypeDeserializer('tests/test_fs_meta.yaml')
        fs_projector = Projector(driver, projection_configuration.root_projection_uri,
                                 projection_configuration.prototype_tree)
        created_projections = [n.get_path() for n in fs_projector.projection_tree.get_tree_nodes()]

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(10, len(created_projections),
                         msg='FS projector created {0} projections.'.format(len(created_projections)))

        # Check root dir projection creation
        self.assertIn('/test_dir',
                      created_projections,
                      msg='Checking creation of root projection')

        # Checking if FASTA files not present in config are not projected
        for i in range(1,5):
            fasta_proj_path = '/test_dir/fasta_file_{0}.fasta'.format(i)
            self.assertNotIn(fasta_proj_path, created_projections,
                          msg='Checking if {0} not in projections'.format(fasta_proj_path))

        # Check BAM files projections creation
        for i in range(1,4):
            bam_proj_path = '/test_dir/bam_file_{0}.bam'.format(i)
            self.assertIn(bam_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(bam_proj_path))

        # Check BAM files metadata projections creation
        for i in range(1,4):
            bam_meta_proj_path = '/test_dir/bam_file_{0}_metadata.json'.format(i)
            self.assertIn(bam_meta_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(bam_meta_proj_path))

        # Check BAM files projections filtration according to metadata
        for i in range(4,6):
            bam_proj_path = '/test_dir/bam_file_{0}.bam'.format(i)
            self.assertNotIn(bam_proj_path, created_projections,
                          msg='Checking creation of {0} projection'.format(bam_proj_path))
