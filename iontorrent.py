#!/usr/bin/env python3
__author__ = 'abragin'

import logging
import io
import json
import os
import stat
import sys
import time
import urllib.request
from urllib.parse import urljoin

from projections import Projection, ProjectionManager
from filesystem import ProjectionFilesystem
from fuse import FUSE


logger = logging.getLogger('iontorrent_projection')

class IonTorrentProjection(ProjectionManager):

    def __init__(self, host, user, password):
        logger.info('Creating Ion Torrent projection for host: %s', host)
        self.host_url = 'http://{}'.format(host)
        self.api_url = 'http://{}/rundb/api/v1/'.format(host)
        self.files_url = urljoin(self.host_url, '/auth/output/Home/')
        self.authenticate(user, password)

        # TODO: switch to tree-like structure instead of manual path parsing
        self.projections = {}
        self.create_projections()

    def authenticate(self, user, password):
        password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        password_manager.add_password(None, self.host_url, user, password)

        handler = urllib.request.HTTPBasicAuthHandler(password_manager)

        # create "opener" (OpenerDirector instance)
        opener = urllib.request.build_opener(handler)

        # use the opener to fetch a URL
        opener.open(self.api_url)

        # Install the opener.
        # Now all calls to urllib.request.urlopen use our opener.
        urllib.request.install_opener(opener)

    def create_projections(self):
        """
        STUB implementation that requests predefined number of experiments.

        :return: list of projections
        """
        # Select last five experiments that were finished (with 'run' status)
        with urllib.request.urlopen(self.api_url+'experiment?status=run&limit=1&order_by=-id') as f:
            experiments = json.loads(f.readall().decode('utf-8'))
        logger.info('Got experiments data: %s', len(experiments['objects']))

        projections = []

        for o in experiments['objects']:
            # Create experiment directory projection
            logger.debug('Got experiments with id: %s, name: %s', o['id'], o['displayName'])

            projection = Projection('/' + o['displayName'], self.host_url+o['resource_uri'])
            projection.type = stat.S_IFDIR

            logger.debug('Created experiment projection: %s', projection)
            projections.append(projection)

            # Create experiment metadata file projection
            exp_meta_projection = Projection('/' + os.path.join(o['displayName'], 'metadata.json'), urljoin(self.api_url, o['resource_uri']))
            logger.debug('Created experiment metadata projection: %s', exp_meta_projection)
            exp_meta_projection.type = stat.S_IFREG
            projections.append(exp_meta_projection)

            # Create experiment plan metadata projection
            plannedexp_metadata_projection = Projection('/' + os.path.join(o['displayName'], 'plannedexperiment.json'), urljoin(self.api_url, o['plan']))
            logger.debug('Created planned experiment metadata projection: %s', plannedexp_metadata_projection)
            projections.append(plannedexp_metadata_projection)

            if o['eas_set']:
                barcoded_samples = o['eas_set'][0]['barcodedSamples']
                logger.debug('Barcoded samples found: %s', barcoded_samples)

                barcodes = {}
                for b in barcoded_samples:
                    barcodes[b] = barcoded_samples[b]['barcodes']
                logger.debug('Barcodes: %s', barcodes)

            # Get sample barcodes data
            with urllib.request.urlopen(self.host_url+o['plan']) as p:
                ex_plan = json.loads(p.readall().decode('utf-8'))
                sample_barcodes = ex_plan['barcodedSamples']

            # Create experiment results directory projections
            for r in o['results']:
                with urllib.request.urlopen(self.host_url + r) as f:
                    results = json.loads(f.readall().decode('utf-8'))

                    path_to_files = os.path.basename(results['filesystempath'])

                    path_to_results_dir = os.path.join(o['displayName'], path_to_files)
                    results_dir_projection = Projection('/'+path_to_results_dir, urljoin(self.files_url, path_to_files))
                    results_dir_projection.type = stat.S_IFDIR

                    projections.append(results_dir_projection)
                    results_metadata_projection = Projection(os.path.join('/'+path_to_results_dir, path_to_files+'.json'), self.host_url+results['resource_uri'])
                    projections.append(results_metadata_projection)

                # Dict stores each variant calling run with result variant call directory as key
                variant_calls = dict()
                with urllib.request.urlopen(urljoin(self.api_url, 'pluginresult?result={}'.format(results['id']))) as f:
                    plugin_res = json.loads(f.readall().decode('utf-8'))
                    for p in plugin_res['objects']:
                        if 'variantCaller' in p['pluginName'] and not 'VFNA' in p['pluginName']:
                            variant_calls[p['path']] = {'barcodes': results['pluginStore'][p['pluginName']]['barcodes'].keys(),
                                                        'resource_uri': p['resource_uri']}
                            # If there is bed file in "target_bed" field, create it`s projection
                            if 'targets_bed' in results['pluginStore'][p['pluginName']]:
                                # Bed file base name
                                bed_file_name = os.path.basename(results['pluginStore'][p['pluginName']]['targets_bed'])
                                # Setting up bed file URI
                                bed_file_path = os.path.join(path_to_files, 'plugin_out',
                                                        os.path.basename(p['path']), bed_file_name)
                                bed_file_projection = Projection(os.path.join('/'+path_to_results_dir, os.path.basename(bed_file_path)),
                                                                 urljoin(self.files_url, bed_file_path))
                                projections.append(bed_file_projection)

                # Create samples projections
                for s in o['samples']:
                    s_path = os.path.join(path_to_results_dir, s['name'])
                    # Creating sample directory projection
                    s_projection = Projection('/' + s_path, urljoin(self.api_url, s['resource_uri']))
                    s_projection.type = stat.S_IFDIR

                    # Getting sample barcode
                    sample_barcode = sample_barcodes[s['name']]['barcodes'][0]
                    # Adding sample BAM file suffix
                    sample_bam_name = sample_barcode + '_rawlib.bam'
                    # Joining sample uri
                    sample_bam_uri = urljoin(self.files_url, os.path.join(path_to_files, sample_bam_name))
                    # Creating sample BAM file projection
                    s_bam_projection = Projection(os.path.join(s_projection.path, s['name'] + '.bam'), sample_bam_uri)
                    # Creating sample metadata projection
                    s_meta_projection = Projection(os.path.join(s_projection.path, 'metadata.json'), urljoin(self.api_url, s['resource_uri']))

                    for vc_path, item in variant_calls.items():
                        # Joining path to variant call directory for sample
                        vc_dir_path = os.path.join(s_projection.path, os.path.basename(vc_path))
                        # Creating variant call directory
                        vc_dir_projection = Projection(vc_dir_path, item['resource_uri'])
                        vc_dir_projection.type = stat.S_IFDIR
                        projections.append(vc_dir_projection)
                        # Joining variant call settings URI
                        vc_settings_path = os.path.join(path_to_files, 'plugin_out',
                                                        os.path.basename(vc_path),
                                                        'local_parameters.json')
                        # Creating variant call settings projection in variant call dir
                        vc_settings_projection = Projection(os.path.join(vc_dir_path, 'variant_caller_settings.json'),
                                                            urljoin(self.files_url, vc_settings_path))
                        projections.append(vc_settings_projection)
                        # If current sample barcode in variant calling run barcodes
                        if sample_barcode in item['barcodes']:
                            # Setting up base URI to variant calling directory
                            base_vcf_file_projection_path = os.path.join(path_to_files, 'plugin_out',
                                                                    os.path.basename(vc_path),
                                                                    sample_barcode)
                            for variant_file_name in ['TSVC_variants.vcf', 'all.merged.vcf', 'indel_assembly.vcf',
                                                      'indel_variants.vcf', 'small_variants.left.vcf',
                                                      'small_variants.vcf', 'small_variants_filtered.vcf',
                                                      'small_variants.sorted.vcf', 'SNP_variants.vcf']:
                                vc_file_path = os.path.join(base_vcf_file_projection_path, variant_file_name)
                                vc_file_projection = Projection(os.path.join(vc_dir_path, variant_file_name),
                                                                 urljoin(self.files_url, vc_file_path))
                                projections.append(vc_file_projection)

                    logging.debug('Created sample projection: %s', s_projection)
                    projections.append(s_projection)
                    projections.append(s_bam_projection)
                    projections.append(s_meta_projection)

        for p in projections:
            self.projections[p.path] = p

    def is_managing_path(self, path):
        return path in self.projections

    def get_projections(self, path):
        logger.info('Requesting projections for path: %s', path)
        projections = []
        for p in self.projections:
            if p.startswith(path):
                # limit projections to one level only
                logging.debug('Analyzing projection candidate: %s', p)
                suffix = p[len(path):]
                if suffix and suffix[0] == '/':
                    suffix = suffix[1:]
                logger.debug('Path suffix: %s', suffix)
                if suffix and '/' not in suffix:
                    projections.append('/' + suffix)

        logger.debug('Returning projections: %s', projections)
        return projections

    def get_attributes(self, path):
        assert path in self.projections

        projection = self.projections[path]

        now = time.time()
        attributes = dict()

        # Set projection attributes

        # This is implementation specific and should be binded to projector data
        attributes['st_atime'] = now
        # This may be implemented as last projection cashing time is casing is enabled
        attributes['st_mtime'] = now
        # On Unix this is time for metedata modification we can use the same conception
        attributes['st_ctime'] = now
        # If this is projection the size is zero
        attributes['st_size'] = projection.size
        # Set type to link anf grant full access to everyone
        attributes['st_mode'] = (projection.type | 0o0777)
        # Set number of hard links to 0
        attributes['st_nlink'] = 0
        # Set id as inode number.
        attributes['st_ino'] = 1

        return attributes

    def open_resource(self, path):
        uri = self.projections[path].uri

        with urllib.request.urlopen(uri) as f:
            content = f.readall()
        logger.info('Got path content: %s\n', path)

        self.projections[path].size = len(content)

        file_header = 3
        resource_io = io.BytesIO(content)

        return file_header, resource_io


# For smoke testing
def main(mountpoint, data_folder, foreground=True):
    # Specify FUSE mount options as **kwargs here. For value options use value=True form, e.g. nonempty=True
    # For complete list of options see: http://blog.woralelandia.com/2012/07/16/fuse-mount-options/
    projection_filesystem = ProjectionFilesystem(mountpoint, data_folder)
    projection_filesystem.projection_manager = IonTorrentProjection('10.5.20.13', 'ionadmin', 'ionadmin')
    fuse = FUSE(projection_filesystem, mountpoint, foreground=foreground, nonempty=True)
    return fuse

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    logging.config.fileConfig(os.path.join(script_dir, 'logging.cfg'))

    # TODO: replace with argparse
    # Get mount point from args
    if len(sys.argv) != 3:
        print('usage: %s <mountpoint> <data folder>' % sys.argv[0])
        exit(1)

    main(sys.argv[1], sys.argv[2])




