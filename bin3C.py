from mzd.cluster import *
from mzd.contact_map import *
from mzd.exceptions import ApplicationException
from mzd.io_utils import load_object, save_object
from mzd.utils import *
import logging
import sys

__version__ = '0.1'

if __name__ == '__main__':
    import argparse

    def mk_version():
        return 'bin3C v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-V', '--version', help='Show the application version')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--clobber', default=False, action='store_true', help='Clobber existing files')
    global_parser.add_argument('--log', help='Log file path [OUTDIR/bin3C.log]')
    global_parser.add_argument('--max-image', type=int, default=4000, help='Maximum image size for plots [4000]')
    global_parser.add_argument('--min-extent', type=int, default=50000,
                               help='Minimum cluster extent used in output [50000]')
    global_parser.add_argument('--min-reflen', type=int, default=1000,
                               help='Minimum acceptable reference length [1000]')
    global_parser.add_argument('--min-signal', type=int, default=5, help='Minimum acceptable signal [5]')

    parser = argparse.ArgumentParser(description='bin3C: a Hi-C based metagenome deconvolution tool')
    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')

    cmd_mkmap = subparsers.add_parser('mkmap', parents=[global_parser],
                                      description='Create a new contact map from assembly sequences and Hi-C bam file.')
    cmd_cluster = subparsers.add_parser('cluster', parents=[global_parser],
                                        description='Cluster an existing contact map into genome bins.')

    """
    make and save the contact map object
    """
    cmd_mkmap.add_argument('-s', '--seed', default=None, help='Random seed')
    cmd_mkmap.add_argument('--eta', default=False, action='store_true',
                           help='Pre-count bam alignments to provide an ETA')
    cmd_mkmap.add_argument('--bin-size', type=int,
                           help='Size of bins for windows extent maps [disabled]')
    cmd_mkmap.add_argument('--min-insert', type=int,
                           help='Minimum pair separation [None]')
    cmd_mkmap.add_argument('--min-mapq', type=int, default=60,
                           help='Minimum acceptable mapping quality [60]')
    cmd_mkmap.add_argument('--strong', type=int, default=10,
                           help='Accepted alignments must being N matches [10]')
    cmd_mkmap.add_argument('-e', '--enzyme', metavar='NEB_NAME', required=True, action='append',
                           help='Case-sensitive NEB enzyme name. Use multiple times for multiple enzymes')
    cmd_mkmap.add_argument('FASTA', help='Reference fasta sequence')
    cmd_mkmap.add_argument('BAM', help='Input bam file in query order')
    cmd_mkmap.add_argument('OUTDIR', help='Output directory')

    """
    cluster the map and save results
    """
    cmd_cluster.add_argument('-s', '--seed', default=None, help='Random seed')
    cmd_cluster.add_argument('--no-report', default=False, action='store_true',
                             help='Do not generate a cluster report')
    cmd_cluster.add_argument('--no-spades', default=False, action='store_true',
                             help='Assembly was not done using SPAdes')
    cmd_cluster.add_argument('--no-plot', default=False, action='store_true',
                             help='Do not generate a clustered heatmap')
    cmd_cluster.add_argument('--no-fasta', default=False, action='store_true',
                             help='Do not generate cluster FASTA files')
    cmd_cluster.add_argument('--only-large', default=False, action='store_true',
                             help='Only write FASTA for clusters longer than min_extent')
    cmd_cluster.add_argument('--algo', default='infomap', choices=['infomap', 'louvain', 'mcl', 'slm', 'simap'],
                             help='Clustering algorithm to apply [infomap]')
    cmd_cluster.add_argument('--fasta', default=None,
                             help='Alternative source FASTA location from that supplied during mkmap')
    cmd_cluster.add_argument('MAP', help='Contact map')
    cmd_cluster.add_argument('OUTDIR', help='Output directory')

    args = parser.parse_args()

    if args.version:
        print mk_version()
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.clobber)
    except IOError as e:
        print 'Error: {}'.format(e.message)
        sys.exit(1)

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
        log_path = os.path.join(args.OUTDIR, 'bin3C.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))

    try:

        if args.command == 'mkmap':

            # Create a contact map for analysis
            cm = ContactMap(args.BAM,
                            args.enzyme,
                            args.FASTA,
                            args.min_insert,
                            args.min_mapq,
                            min_len=args.min_reflen,
                            min_sig=args.min_signal,
                            min_extent=args.min_extent,
                            strong=args.strong,
                            bin_size=args.bin_size,
                            precount=args.eta)

            if cm.is_empty():
                logger.info('Stopping as the map is empty')
                sys.exit(1)

            logger.info('Saving contact map instance')
            save_object(os.path.join(args.OUTDIR, 'contact_map.p'), cm)

        elif args.command == 'cluster':

            if not args.seed:
                args.seed = make_random_seed()
                logger.info('Generated random seed: {}'.format(args.seed))
            else:
                logger.info("User set random seed: {}".format(args.seed))

            # Load a pre-existing serialized contact map
            logger.info('Loading existing contact map from: {}'.format(args.MAP))
            cm = load_object(args.MAP)
            cm.min_extent = args.min_extent
            # update the mask if the user has changed the thresholds
            if args.min_signal != cm.min_sig or args.min_reflen != cm.min_len:
                # pedantically set these and pass to method just in-case of logic oversight
                cm.min_len = args.min_reflen
                cm.min_sig = args.min_signal
                cm.set_primary_acceptance_mask(min_sig=args.min_signal, min_len=args.min_reflen, update=True)

            # cluster the entire map
            clustering = cluster_map(cm, method='infomap', seed=args.seed, work_dir=args.OUTDIR)
            # generate report per cluster
            cluster_report(cm, clustering, is_spades=not args.no_spades)
            # write MCL clustering file
            write_mcl(cm, os.path.join(args.OUTDIR, 'clustering.mcl'), clustering)
            # serialize full clustering object
            save_object(os.path.join(args.OUTDIR, 'clustering.p'), clustering)

            if not args.no_report:
                # write a tabular report
                write_report(os.path.join(args.OUTDIR, 'cluster_report.csv'), clustering)

            if not args.no_fasta:
                # write per-cluster fasta files, also separate ordered fasta if an ordering exists
                write_fasta(cm, args.OUTDIR, clustering, source_fasta=args.fasta, clobber=True,
                            only_large=args.only_large)

            if not args.no_plot:
                # the entire clustering
                plot_clusters(cm, os.path.join(args.OUTDIR, 'cluster_plot.png'), clustering,
                              max_image_size=args.max_image, ordered_only=False, simple=False, permute=True)

    except ApplicationException as ex:
        import sys
        logger.error(ex.message)
        sys.exit(1)
