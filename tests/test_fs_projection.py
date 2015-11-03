import os
import time
import subprocess
import logging
import logging.config
from unittest import TestCase, skip
from projections import PrototypeDeserializer, Projection
from fs_projection import FSDriver, FSProjector

MOUNT_POINT = 'tests/mnt'
DATA_FOLDER = 'tests/data'
CONFIG_PATH = 'tests/test_fs_proj_config.yaml'

# Import logging configuration from the file provided
logging.config.fileConfig('logging.cfg')
logger = logging.getLogger('sra_test')


class TestFSProjection(TestCase):

    def setUp(self):
        driver = FSDriver()
        projection_configuration = PrototypeDeserializer(CONFIG_PATH)
        root_projection = Projection('/', projection_configuration.root_projection_uri)
        self.fs_projector = FSProjector(driver, root_projection, projection_configuration.prototype_tree)

    def test_create_projections(self):
        """
        Tests if filesystem projector creates projections
        """

        created_projections = self.fs_projector.projections

        # Test if number of created projections equals to expected number of projections
        self.assertEqual(37, len(created_projections),
                         msg='Checking if FS projector created 37 projections, current number: {0}'.format(len(created_projections)))

        # Check root dir projection creation
        self.assertIn('/test_dir',
                      created_projections,
                      msg='Checking creation of root projection')

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
        # Unmount any previous tests
        subprocess.Popen(['fusermount', '-u', MOUNT_POINT])

        fs_proj = subprocess.Popen(['./fs_projection.py',
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

