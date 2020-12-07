from proxigenomics_toolkit.contact_map import *
from proxigenomics_toolkit.exceptions import ApplicationException
from proxigenomics_toolkit.io_utils import load_object, save_object
from proxigenomics_toolkit.misc_utils import make_random_seed
from bin3C._version import version_stamp

import argparse
import logging
import sys
from pipes import quote


def reconstruct_cmdline():
    """
    Reconstruct what could have been the command line, where arguments are properly escaped
    and quoted.
    :return: string representation of command line options
    """
    return ' '.join(map(quote, sys.argv))


def or_default(v, default):
    if v is not None:
        assert isinstance(v, type(default)), \
            'supplied value [{}] is not of the correct type [{}]'.format(v, type(default))
        return v
    return default


def required_length(n_min):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if len(values) < n_min:
                msg = 'argument "{f}" requires at least {nmin} arguments'.format(f=self.dest, nmin=n_min)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


def main():

    runtime_defaults = {
        'min_reflen': 1000,
        'min_signal': 2,
        'max_image': 4000,
        'min_extent': 50000,
        'min_mapq': 60,
        'max_edist': 2,
        'min_alen': 25,
    }

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--clobber', default=False, action='store_true', help='Clobber existing files')
    global_parser.add_argument('--log', help='Log file path [OUTDIR/bin3C.log]')

    parser_top = argparse.ArgumentParser(description='bin3C: a Hi-C based metagenome deconvolution tool',
                                         add_help=False)
    parser_top.add_argument('-V', '--version', default=False, action='store_true', help='Version')

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')
    subparsers.required = False
    cmd_mkmap = subparsers.add_parser('mkmap', parents=[global_parser],
                                      description='Create a new contact map from assembly sequences and Hi-C bam file.')
    cmd_combine = subparsers.add_parser('combine', parents=[global_parser],
                                        description='Additively combine contact maps derived from the same assembly.')
    cmd_cluster = subparsers.add_parser('cluster', parents=[global_parser],
                                        description='Cluster an existing contact map into genome bins.')
    cmd_extract = subparsers.add_parser('extract', parents=[global_parser],
                                        description='Extract a representation of single cluster')

    """
    make and save the contact map object
    """
    cmd_mkmap.add_argument('--threads', metavar='INT', type=int, default=1,
                           help='Number of IO threads used when accessing BAM files [1]')
    cmd_mkmap.add_argument('--bin-size', metavar='NBASES', type=int,
                           help='Size of bins for windows extent maps [None]')
    cmd_mkmap.add_argument('--keep-duplicates', default=False, action='store_true',
                           help='Do not remove duplicate pair mappings')
    cmd_mkmap.add_argument('--eta', default=False, action='store_true',
                           help='Count bam alignments to provide an estimated processing time')
    cmd_mkmap.add_argument('--min-extent', metavar='NBASES', type=int,
                           help='Minimum cluster extent used in output [50000]')
    cmd_mkmap.add_argument('--min-reflen', metavar='NBASES', type=int,
                           help='Minimum acceptable reference length [1000]')
    cmd_mkmap.add_argument('--min-signal', metavar='COUNTS', type=int,
                           help='Minimum acceptable signal [2]')
    cmd_mkmap.add_argument('--min-insert', metavar='NBASES', type=int,
                           help='Minimum pair separation [None]')
    cmd_mkmap.add_argument('--min-mapq', metavar='INT', type=int,
                           help='Minimum acceptable mapping quality [60]')
    cmd_mkmap.add_argument('--max-edist', metavar='INT', type=int,
                           help='Maximum acceptable edit distance [2]')
    cmd_mkmap.add_argument('--min-alen', metavar='NBASES', type=int,
                           help='Minimum acceptable alignment length [25]')
    cmd_mkmap.add_argument('--tip-size', metavar='NBASES', type=int, default=None,
                           help='The size of the region used when tracking only the ends '
                                'of contigs [Experimental] [None]')
    cmd_mkmap.add_argument('-e', '--enzyme', metavar='NEB_NAME', required=True, action='append',
                           help='Case-sensitive NEB enzyme name. Use multiple times for multiple enzymes')
    cmd_mkmap.add_argument('FASTA', help='Reference fasta sequence')
    cmd_mkmap.add_argument('BAM', help='Input bam file in query order')
    cmd_mkmap.add_argument('OUTDIR', help='Output directory')

    """
    combine contact maps
    """
    cmd_combine.add_argument('OUTDIR', help='Output directory')
    cmd_combine.add_argument('MAP', nargs='+', help='Contact maps to combine', action=required_length(2))

    """
    cluster the map and save results
    """
    cmd_cluster.add_argument('-s', '--seed', metavar='INT', default=None, help='Random seed')
    cmd_cluster.add_argument('--min-extent', metavar='NBASES', type=int,
                             help='Minimum cluster extent used in output [50000]')
    cmd_cluster.add_argument('--min-reflen', metavar='NBASES', type=int,
                             help='Minimum acceptable reference length [1000]')
    cmd_cluster.add_argument('--min-signal', metavar='COUNTS', type=int,
                             help='Minimum acceptable signal [2]')
    cmd_cluster.add_argument('--max-image', metavar='PIXELS', type=int, help='Maximum image size for plots [4000]')
    cmd_cluster.add_argument('--no-report', default=False, action='store_true',
                             help='Do not generate a cluster report')
    cmd_cluster.add_argument('--assembler', choices=['generic', 'spades', 'megahit'], default='generic',
                             help='Assembly software used to create contigs')
    cmd_cluster.add_argument('--no-plot', default=False, action='store_true',
                             help='Do not generate a clustered heatmap')
    cmd_cluster.add_argument('--plot-format', default='png', choices=['png', 'pdf'],
                             help='File format when writing contact map plot [png]')
    cmd_cluster.add_argument('--norm-method', default='sites', choices=['sites', 'gothic-es', 'gothic-sig'],
                             help='Contact map normalisation method [sites]')
    cmd_cluster.add_argument('--no-fasta', default=False, action='store_true',
                             help='Do not generate cluster FASTA files')
    cmd_cluster.add_argument('--only-large', default=False, action='store_true',
                             help='Only write FASTA for clusters longer than min_extent')
    cmd_cluster.add_argument('--coverage', metavar='PATH', default=None,
                             help='Per-sequence depth of coverage data format: "seq_id,value"')
    # cmd_cluster.add_argument('--algo', default='infomap', choices=['infomap', 'louvain', 'mcl', 'slm', 'simap'],
    #                          help='Clustering algorithm to apply [infomap]')
    cmd_cluster.add_argument('--fasta', metavar='PATH', default=None,
                             help='Alternative location of source FASTA from that supplied during mkmap')
    cmd_cluster.add_argument('--n-iter', '-N', metavar="INT", default=None, type=int,
                             help='Number of iterations for clustering optimisation [10]')
    cmd_cluster.add_argument('MAP', help='bin3C contact map')
    cmd_cluster.add_argument('OUTDIR', help='Output directory')

    """
    extract a single cluster
    """

    cmd_extract.add_argument('--threads', metavar='INT', type=int, default=1,
                             help='Number of IO threads for accessing BAM file [1]')
    cmd_extract.add_argument('--max-image', metavar='PIXELS', type=int, help='Maximum image size for plots [4000]')
    cmd_extract.add_argument('--use-extent', default=False, action='store_true',
                             help='For plots use extent map rather than sequence map if available')
    cmd_extract.add_argument('--show-sequences', default=False, action='store_true',
                             help='For plots grid lines and labels mark individual '
                                  'sequences rather than whole clusters')
    cmd_extract.add_argument('-b', '--bam', help='Alternative location of source BAM file')
    cmd_extract.add_argument('--plot-format', default='png', choices=['png', 'pdf'],
                             help='File format when writing contact map plot [png]')
    cmd_extract.add_argument('--norm-method', default='sites', choices=['sites', 'gothic-es', 'gothic-sig'],
                             help='Contact map normalisation method [sites]')
    cmd_extract.add_argument('-f', '--format', choices=['graph', 'plot', 'bam'], default='plot',
                             help='Select output format [plot]')
    cmd_extract.add_argument('MAP', help='bin3C contact map')
    cmd_extract.add_argument('CLUSTERING', help='bin3C clustering object')
    cmd_extract.add_argument('OUTDIR', help='Output directory')
    cmd_extract.add_argument('CLUSTER_ID', nargs='*', type=int, help='1-based Cluster number (eg. 1,2,..,99)')

    args, extras = parser_top.parse_known_args()

    if args.version:
        print version_stamp()
        sys.exit(0)

    if len(extras) == 0:
        parser.print_usage()
        sys.exit(0)

    parser.parse_args(extras, namespace=args)

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
    logger.debug(version_stamp(False))
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(reconstruct_cmdline()))

    try:

        if args.command == 'mkmap':

            if args.tip_size is not None:
                logger.warn('[Experimental] Specifying tip-size enables independent tracking of the ends of contigs.')
                if args.tip_size < 5000:
                    logger.warn('[Experimental] It is recommended to use tip sizes no smaller than 5kbp')
                if args.tip_size > args.min_reflen:
                    msg = 'min-reflen cannot be smaller than the tip-size'
                    logger.error('[Experimental] {}'.format(msg))
                    raise ApplicationException(msg)

            # Create a contact map for analysis
            cm = ContactMap(args.BAM,
                            args.enzyme,
                            args.FASTA,
                            args.min_insert,
                            or_default(args.min_mapq, runtime_defaults['min_mapq']),
                            min_len=or_default(args.min_reflen, runtime_defaults['min_reflen']),
                            min_sig=or_default(args.min_signal, runtime_defaults['min_signal']),
                            min_extent=or_default(args.min_extent, runtime_defaults['min_extent']),
                            max_edist=or_default(args.max_edist, runtime_defaults['max_edist']),
                            min_alen=or_default(args.min_alen, runtime_defaults['min_alen']),
                            bin_size=args.bin_size,
                            tip_size=args.tip_size,
                            no_duplicates=not args.keep_duplicates,
                            precount=args.eta,
                            threads=args.threads)

            if cm.is_empty():
                logger.info('Stopping as the map is empty')
                sys.exit(1)

            logger.info('Saving contact map instance')
            save_object(os.path.join(args.OUTDIR, 'contact_map.p'), cm)

        elif args.command == 'combine':
            cm_summed = None
            for cm_file in args.MAP:
                logger.info('Loading contact map: {}'.format(cm_file))
                if cm_summed is None:
                    cm_summed = load_object(cm_file)
                else:
                    cm_summed.append_map(load_object(cm_file))
            logger.info('Saving combined map')
            save_object(os.path.join(args.OUTDIR, 'combined_map.p'), cm_summed)

        elif args.command == 'cluster':

            if not args.seed:
                args.seed = make_random_seed()
                logger.info('Generated random seed: {}'.format(args.seed))
            else:
                logger.info("User set random seed: {}".format(args.seed))

            # Load a pre-existing serialized contact map
            logger.info('Loading existing contact map from: {}'.format(args.MAP))
            cm = load_object(args.MAP)

            if args.min_extent is not None:
                cm.min_extent = args.min_extent

            # in cases where a user supplies a value, we will need to redo the acceptance mask
            # otherwise the mask will have been done with default values.
            remask = False
            if args.min_signal is not None:
                cm.min_sig = args.min_signal
                remask = True
            if args.min_reflen is not None:
                cm.min_len = args.min_reflen
                remask = True

            if remask:
                cm.set_primary_acceptance_mask(min_sig=cm.min_sig, min_len=cm.min_len, update=True)

            # cluster the entire map
            clustering = cluster_map(cm, method='infomap', seed=args.seed, work_dir=args.OUTDIR,
                                     n_iter=args.n_iter, norm_method=args.norm_method)
            if not args.no_report:
                # generate report per cluster
                cluster_report(cm, clustering, assembler=args.assembler, source_fasta=args.fasta,
                               coverage_file=args.coverage)
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
                plot_clusters(cm, os.path.join(args.OUTDIR, 'cluster_plot.{}'.format(args.plot_format)), clustering,
                              max_image_size=or_default(args.max_image, runtime_defaults['max_image']),
                              ordered_only=False, simple=False, permute=True)

        elif args.command == 'extract':

            logger.info('Loading contact map from: {}'.format(args.MAP))
            cm = load_object(args.MAP)

            logger.info('Loading clustering solution from: {}'.format(args.CLUSTERING))
            clustering = load_object(args.CLUSTERING)

            # Convert public string ids to internal 0-based integer ids
            if not args.CLUSTER_ID:
                cluster_ids = None
                logger.info('Extracting all clusters')
            else:
                cluster_ids = np.asarray(args.CLUSTER_ID, dtype=np.int) - 1
                logger.info('Extracting {} clusters'.format(len(cluster_ids)))

            if args.format in ['plot', 'graph']:
                # ensure that selected clusters are not masked
                cm.min_sig = 0
                cm.min_len = 0
                cm.min_extent = 0
                cm.set_primary_acceptance_mask(min_sig=cm.min_sig, min_len=cm.min_len, update=True)
                cm.prepare_seq_map(norm=True, bisto=True, norm_method=args.norm_method)

            if args.format == 'plot':

                if args.use_extent and cm.extent_map is None:
                    logger.error('An extent map was not generated when creating the specified contact map')

                plot_clusters(cm, os.path.join(args.OUTDIR, 'extracted.{}'.format(args.plot_format)), clustering,
                              max_image_size=or_default(args.max_image, runtime_defaults['max_image']),
                              ordered_only=False, permute=True, simple=not args.use_extent,
                              cl_list=cluster_ids, norm_method=args.norm_method, show_sequences=args.show_sequences)

            elif args.format == 'graph':

                g = to_graph(cm, norm=True, bisto=True, node_id_type='external', scale=False,
                             clustering=clustering, cl_list=cluster_ids, norm_method=args.norm_method)

                nx.write_graphml(g, os.path.join(args.OUTDIR, 'extracted.graphml'))

            elif args.format == 'bam':

                out_file, n_refs, n_pairs = extract_bam(cm, clustering, args.OUTDIR, cluster_ids, clobber=args.clobber,
                                                        threads=args.threads, bam_file=args.bam,
                                                        version=version_stamp(False), cmdline=reconstruct_cmdline())
                logger.info('Output BAM {} contains {:,} references and {:,} pairs'.format(out_file, n_refs, n_pairs))

            else:
                raise ApplicationException('Unknown format option {}'.format(args.format))

    except ApplicationException as ex:
        logger.error(ex.message)
        sys.exit(1)
