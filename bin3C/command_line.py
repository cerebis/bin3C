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


def required_length(n_min):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if len(values) < n_min:
                msg = 'argument "{f}" requires at least {nmin} arguments'.format(f=self.dest, nmin=n_min)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


def main():

    _defaults = {
        'min_reflen': 1000,
        'min_signal': 2,
        'max_image': 4000,
        'min_extent': 5000,
        'min_insert': None,
        'min_mapq': 60,
        'max_edist': 4,
        'min_alen': 25,
        'threads': 1,
        'bin_size': None,
        'tip_size': None,
        'n-iter': 10,
        'norm-method': 'sites',
    }

    # options shared by all commands
    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-d', '--debug', default=False, action='store_true',
                               help='Fall into debug mode on exception')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--clobber', default=False, action='store_true', help='Clobber existing files')
    global_parser.add_argument('--log', help='Log file path (default: OUTDIR/bin3C.log)')

    # options shared by primary analysis commands
    analysis_parser = argparse.ArgumentParser(add_help=False)
    analysis_parser.add_argument('--min-extent', metavar='NBASES', type=int, default=_defaults['min_extent'],
                                 help='Minimum cluster extent used in output (default: %(default)s)')
    analysis_parser.add_argument('--min-reflen', metavar='NBASES', type=int, default=_defaults['min_reflen'],
                                 help='Minimum acceptable reference length (default: %(default)s)')
    analysis_parser.add_argument('--min-signal', metavar='COUNTS', type=int, default=_defaults['min_signal'],
                                 help='Minimum acceptable signal (default: %(default)s)')

    parser = argparse.ArgumentParser(description='bin3C: a tool for Hi-C based metagenome-assembled genome binning',
                                     add_help=True)
    parser.add_argument('-V', '--version', default=False, action='store_true', help='Version')

    # sub-commands fall beneath this parser
    command_parsers = parser.add_subparsers(title='Valid subcommands', dest='command',
                                            help='Specify a subcommand for further options')
    command_parsers.required = False

    """
    make and save the contact map object
    """
    cmd_mkmap = command_parsers.add_parser('mkmap', parents=[global_parser, analysis_parser],
                                           description='Create a contact map from assembly '
                                                       'sequences and Hi-C bam file.')
    cmd_mkmap.add_argument('--threads', metavar='INT', type=int, default=_defaults['threads'],
                           help='Number of IO threads used when accessing BAM files (default: %(default)s)')
    cmd_mkmap.add_argument('--bin-size', metavar='NBASES', type=int, default=None,
                           help='Size of bins for windows extent maps (default: %(default)s)')
    cmd_mkmap.add_argument('--keep-duplicates', default=False, action='store_true',
                           help='Do not remove duplicate pair mappings')
    cmd_mkmap.add_argument('--eta', default=False, action='store_true',
                           help='Count bam alignments to provide an estimated processing time')
    cmd_mkmap.add_argument('--min-insert', metavar='NBASES', type=int, default=_defaults['min_insert'],
                           help='Minimum pair separation (default: %(default)s)')
    cmd_mkmap.add_argument('--min-mapq', metavar='INT', type=int, default=_defaults['min_mapq'],
                           help='Minimum acceptable mapping quality (default: %(default)s)')
    cmd_mkmap.add_argument('--max-edist', metavar='INT', type=int, default=_defaults['max_edist'],
                           help='Maximum acceptable edit distance (default: %(default)s)')
    cmd_mkmap.add_argument('--min-alen', metavar='NBASES', type=int, default=_defaults['min_alen'],
                           help='Minimum acceptable alignment length (default: %(default)s)')
    cmd_mkmap.add_argument('--tip-size', metavar='NBASES', type=int, default=_defaults['tip_size'],
                           help='The size of the region used when tracking only the ends '
                                'of contigs [Experimental] (default: %(default)s)')
    digestion_group = cmd_mkmap.add_mutually_exclusive_group(required=True)
    digestion_group.add_argument('-k', '--library-kit', choices=['phase', 'arima'], default=None, nargs='?',
                                 help='Define enzymes by commercial library kit')
    digestion_group.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append',
                                 help='Case-sensitive NEB enzyme name. Use multiple times for multiple enzymes')
    cmd_mkmap.add_argument('FASTA', help='Reference fasta sequence')
    cmd_mkmap.add_argument('BAM', help='Input bam file in query order')
    cmd_mkmap.add_argument('OUTDIR', help='Output directory')

    """
    cluster the map and save results
    """
    cmd_cluster = command_parsers.add_parser('cluster',
                                             parents=[global_parser, analysis_parser],
                                             description='Cluster an existing contact map into genome bins.')
    cmd_cluster.add_argument('-s', '--seed', metavar='INT', default=None, help='Random seed (default: %(default)s)')
    cmd_cluster.add_argument('--max-image', metavar='PIXELS', type=int, default=_defaults['max_image'],
                             help='Maximum image size for plots (default: %(default)s)')
    cmd_cluster.add_argument('--no-report', default=False, action='store_true',
                             help='Do not generate a cluster report')
    cmd_cluster.add_argument('--assembler', choices=['generic', 'spades', 'megahit'], default='generic',
                             help='Assembly software used to create contigs (default: %(default)s)')
    cmd_cluster.add_argument('--no-plot', default=False, action='store_true',
                             help='Do not generate a clustered heatmap')
    cmd_cluster.add_argument('--plot-format', default='png', choices=['png', 'pdf'],
                             help='File format when writing contact map plot (default: %(default)s)')
    cmd_cluster.add_argument('--norm-method', default=_defaults['norm-method'],
                             choices=['sites', 'gothic-effect', 'gothic-binomial', 'gothic-poisson'],
                             help='Contact map normalisation method (default: %(default)s)')
    cmd_cluster.add_argument('--no-fasta', default=False, action='store_true',
                             help='Do not generate cluster FASTA files')
    cmd_cluster.add_argument('--only-large', default=False, action='store_true',
                             help='Only write FASTA for clusters longer than min_extent')
    cmd_cluster.add_argument('--coverage', metavar='PATH', default=None,
                             help='Per-sequence depth of coverage data format: "seq_id,value" (default: %(default)s)')
    # cmd_cluster.add_argument('--algo', default='infomap', choices=['infomap', 'louvain', 'mcl', 'slm', 'simap'],
    #                          help='Clustering algorithm to apply [infomap]')
    cmd_cluster.add_argument('--fasta', metavar='PATH', default=None,
                             help='Alternative location of source FASTA from that supplied during mkmap')
    cmd_cluster.add_argument('--n-iter', '-N', metavar="INT", type=int, default=_defaults['n-iter'],
                             help='Number of iterations for clustering optimisation (default: %(default)s)')
    cmd_cluster.add_argument('--exclude-from', metavar='FILE', default=None,
                             help='File of sequence ids (one-per-line) to '
                                  'exclude from clustering (default: %(default)s)')
    cmd_cluster.add_argument('MAP', help='bin3C contact map')
    cmd_cluster.add_argument('OUTDIR', help='Output directory')

    """
    combine contact maps
    """
    cmd_combine = command_parsers.add_parser('combine', parents=[global_parser],
                                             description='Additively combine contact maps '
                                                         'derived from the same assembly.')
    cmd_combine.add_argument('OUTDIR', help='Output directory')
    cmd_combine.add_argument('MAP', nargs='+', help='Contact maps to combine', action=required_length(2))

    """
    extract a single cluster
    """
    cmd_extract = command_parsers.add_parser('extract',
                                             parents=[global_parser],
                                             description='Extract a representation of single cluster')
    cmd_extract.add_argument('--threads', metavar='INT', type=int, default=_defaults['threads'],
                             help='Number of IO threads for accessing BAM file (default: %(default)s)')
    cmd_extract.add_argument('--max-image', metavar='PIXELS', type=int, default=_defaults['max_image'],
                             help='Maximum image size for plots (default: %(default)s)')
    cmd_extract.add_argument('--use-extent', default=False, action='store_true',
                             help='For plots use extent map rather than sequence map if available')
    cmd_extract.add_argument('--show-sequences', default=False, action='store_true',
                             help='For plots grid lines and labels mark individual '
                                  'sequences rather than whole clusters')
    cmd_extract.add_argument('-b', '--bam', help='Alternative location of source BAM file')
    cmd_extract.add_argument('--plot-format', default='png', choices=['png', 'pdf'],
                             help='File format when writing contact map plot (default: %(default)s)')
    cmd_extract.add_argument('--norm-method', default=_defaults['norm-method'],
                             choices=['sites', 'gothic-effect', 'gothic-binomial', 'gothic-poisson'],
                             help='Contact map normalisation method (default: %(default)s)')
    cmd_extract.add_argument('-f', '--format', choices=['graph', 'plot', 'bam'], default='plot',
                             help='Select output format (default: %(default)s)')
    cmd_extract.add_argument('MAP', help='bin3C contact map')
    cmd_extract.add_argument('CLUSTERING', help='bin3C clustering object')
    cmd_extract.add_argument('OUTDIR', help='Output directory')
    cmd_extract.add_argument('CLUSTER_ID', nargs='*', type=int, help='1-based Cluster number (eg. 1,2,..,99)')

    args = parser.parse_args()

    if args.version:
        print(version_stamp())
        sys.exit(0)

    if args.command is None:
        parser.print_usage()
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.clobber)
    except IOError as ex:
        print('Error: {}'.format(ex))
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
                logger.warning('[Experimental] Specifying tip-size enables '
                               'independent tracking of the ends of contigs.')
                if args.tip_size < 5000:
                    logger.warning('[Experimental] It is recommended to use tip sizes no smaller than 5kbp')
                if args.tip_size > args.min_reflen:
                    msg = 'min-reflen cannot be smaller than the tip-size'
                    logger.error('[Experimental] {}'.format(msg))
                    raise ApplicationException(msg)

            # check if the user has employed the library kit option to declare enzymes
            if args.command in {'kmer', 'bam'} and args.library_kit is not None:
                # commercial kit definitions
                kit_choices = {'phase': ['Sau3AI', 'MluCI'],
                               'arima': ['DpnII', 'HinfI']}
                args.enzyme = kit_choices[args.library_kit]
                logger.info('Library kit {} declares enzymes {}'.format(args.library_kit, args.enzyme))

            # Create a contact map for analysis
            cm = ContactMap(args.BAM,
                            args.enzyme,
                            args.FASTA,
                            args.min_insert,
                            min_mapq=args.min_mapq,
                            min_len=args.min_reflen,
                            min_sig=args.min_signal,
                            min_extent=args.min_extent,
                            max_edist=args.max_edist,
                            min_alen=args.min_alen,
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

            exclude_names = []
            if args.exclude_from:
                logger.info('Reading excluded ids from {}'.format(args.exclude_from))
                for _nm in open(args.exclude_from, 'r'):
                    _nm = _nm.strip()
                    if not _nm or _nm.startswith('#'):
                        continue
                    exclude_names.append(_nm)

            # cluster the entire map
            clustering = cluster_map(cm,
                                     method='infomap',
                                     seed=args.seed,
                                     work_dir=args.OUTDIR,
                                     n_iter=args.n_iter,
                                     norm_method=args.norm_method,
                                     exclude_names=exclude_names)
            if not args.no_report:
                # generate report per cluster
                cluster_report(cm,
                               clustering,
                               assembler=args.assembler,
                               source_fasta=args.fasta,
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
                write_fasta(cm,
                            args.OUTDIR,
                            clustering,
                            source_fasta=args.fasta,
                            clobber=True,
                            only_large=args.only_large)

            if not args.no_plot:
                # the entire clustering
                plot_clusters(cm,
                              os.path.join(args.OUTDIR, 'cluster_plot.{}'.format(args.plot_format)),
                              clustering,
                              max_image_size=args.max_image,
                              ordered_only=False,
                              simple=False,
                              permute=True)

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

                plot_clusters(cm,
                              os.path.join(args.OUTDIR, 'extracted.{}'.format(args.plot_format)),
                              clustering,
                              max_image_size=args.max_image,
                              ordered_only=False,
                              permute=True,
                              simple=not args.use_extent,
                              cl_list=cluster_ids,
                              norm_method=args.norm_method,
                              show_sequences=args.show_sequences)

            elif args.format == 'graph':

                g = to_graph(cm,
                             norm=True,
                             bisto=True,
                             node_id_type='external',
                             scale=False,
                             clustering=clustering,
                             cl_list=cluster_ids,
                             norm_method=args.norm_method)

                nx.write_graphml(g, os.path.join(args.OUTDIR, 'extracted.graphml'))

            elif args.format == 'bam':

                out_file, n_refs, n_pairs = extract_bam(cm,
                                                        clustering,
                                                        args.OUTDIR,
                                                        cluster_ids,
                                                        clobber=args.clobber,
                                                        threads=args.threads,
                                                        bam_file=args.bam,
                                                        version=version_stamp(False),
                                                        cmdline=reconstruct_cmdline())

                logger.info('Output BAM {} contains {:,} references and {:,} pairs'.format(out_file, n_refs, n_pairs))

            else:
                raise ApplicationException('Unknown format option {}'.format(args.format))

    except ApplicationException as ex:
        logger.error(ex)
        sys.exit(1)
