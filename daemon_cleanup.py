import Pyro4

LOCK_FILE = '/var/lock/projections'


def perform_daemon_cleanup():
    """
    This function performs daemon cleanup termination Projectors subprocessess and flushing daemon`s connection to
    database.
    :return: None
    """
    with open(LOCK_FILE, 'r') as f:
        uri = f.readline()
        projections_daemon = Pyro4.Proxy(uri)

        projections_daemon.stop_daemon()


if __name__ == '__main__':
    perform_daemon_cleanup()
