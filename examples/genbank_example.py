#!/usr/bin/env python3

from projections_daemon import ProjectionsDaemon



def main(cfg_path, mountpoint, data_folder, foreground=True):
    daemon = ProjectionsDaemon()

    daemon.project('genbak_projection', 'mount', 'genbank')


if __name__ == '__main__':
    main()
