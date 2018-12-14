from mzd.splitters import simple_splitter
from mzd.exceptions import ApplicationException
import logging
import sys
import os

__version__ = '0.1'

if __name__ == '__main__':
    import argparse

    def mk_version():
        return 'split v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    def ifelse(arg, default):
        if arg is None:
            return default
        else:
            return arg


    parser = argparse.ArgumentParser(description='Split references prior to mapping Hi-C reads.')

    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    parser.add_argument('--clobber', default=False, action='store_true', help='Clobber existing files')
    parser.add_argument('--log', help='Log file path [split.log]')
    parser.add_argument('-s', '--size', type=int, default=10000,
                           help='The target size in bp for fragments in bp [10000]')
    parser.add_argument('FASTA', help='Input reference fasta sequence')
    parser.add_argument('OUTFILE', help='Output split reference fasta', nargs='?')
    args = parser.parse_args()

    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    if args.log is not None:
        log_path = args.log
    else:
        log_path = os.path.join('split.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    try:

        simple_splitter(args.size, args.FASTA, out_file=args.OUTFILE, clobber=args.clobber)

    except ApplicationException as ex:
        import sys
        logger.error(ex.message)
        sys.exit(1)
