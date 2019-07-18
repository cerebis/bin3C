from . import louvain
from .utils import app_path, make_dir
from .contact_map import SeqOrder
from .seq_utils import IndexedFasta
from .exceptions import *
from typing import Optional
import Bio.SeqIO as SeqIO
import Bio.SeqUtils as SeqUtils
import contextlib
import itertools
import logging
import networkx as nx
import numpy as np
import os
import pandas
import scipy.sparse as sp
import subprocess
import tqdm

logger = logging.getLogger(__name__)


def add_cluster_names(clustering, prefix='CL'):
    """
    Add sequential names beginning from 1 to a clustering in-place.

    A pedantic determination of how many digits are required for the
    largest cluster number is performed, so names will sort conveniently
    in alphanumeric order and align to the eye in output information.

    :param clustering: clustering solution returned from ContactMap.cluster_map
    :param prefix: static prefix of cluster names.
    """
    try:
        num_width = max(1, int(np.ceil(np.log10(max(clustering)+1))))
    except OverflowError:
        num_width = 1

    for cl_id in clustering:
        # names will 1-based
        clustering[cl_id]['name'] = '{0}{1:0{2}d}'.format(prefix, cl_id+1, num_width)


def cluster_map(contact_map, seed, method='infomap', min_len=None, min_sig=None, work_dir='.'):
    """
    Cluster a contact map into groups, as an approximate proxy for "species" bins.

    :param contact_map: an instance of ContactMap to cluster
    :param method: clustering algorithm to employ
    :param seed: a random seed
    :param min_len: override minimum sequence length, otherwise use instance's setting)
    :param min_sig: override minimum off-diagonal signal (in raw counts), otherwise use instance's setting)
    :param work_dir: working directory to which files are written during clustering
    :return: a dictionary detailing the full clustering of the contact map
    """

    def _read_mcl(pathname):
        """
        Read a MCL solution file converting this to a TruthTable.

        :param pathname: mcl file name
        :return: dict of cluster_id to array of seq_ids
        """

        with open(pathname, 'r') as h_in:
            # read the MCL file, which lists all members of a class on a single line
            # the class ids are implicit, therefore we use line number.
            cl_map = {}
            for cl_id, line in enumerate(h_in):
                line = line.rstrip()
                if not line:
                    break
                cl_map[cl_id] = np.array(sorted([int(tok) for tok in line.split()]))
        return cl_map

    def _read_table(pathname, seq_col=0, cl_col=1):
        # type: (str, Optional[int], int) -> dict
        """
        Read cluster solution from a tabular file, one assignment per line. Implicit sequence
        naming is achieved by setting seq_col=None. The reverse (implicit column naming) is
        not currently supported.

        :param pathname: table file name
        :param seq_col: column number of seq_ids
        :param cl_col: column number of cluster ids
        :return: dict of cluster_id to array of seq_ids
        """
        assert seq_col != cl_col, 'sequence and cluster columns must be different'
        with open(pathname, 'r') as h_in:
            cl_map = {}
            n = 0
            for line in h_in:
                line = line.strip()
                if not line:
                    break
                if seq_col is None:
                    cl_id = int(line)
                    seq_id = n
                    n += 1
                else:
                    t = line.split()
                    if len(t) != 2:
                        logger.warning('invalid line encountered when reading cluster table: {}'.format(line))

                    seq_id, cl_id = int(t[seq_col]), int(t[cl_col])
                cl_map.setdefault(cl_id, []).append(seq_id)
            for k in cl_map:
                cl_map[k] = np.array(cl_map[k], dtype=np.int)
            return cl_map

    def _read_tree(pathname):
        """
        Read a tree clustering file as output by Infomap.

        :param pathname: the path to the tree file
        :return: dict of cluster_id to array of seq_ids
        """
        with open(pathname, 'r') as in_h:
            cl_map = {}
            for line in in_h:
                line = line.strip()
                if not line:
                    break
                if line.startswith('#'):
                    continue
                fields = line.split()
                hierarchy = fields[0].split(':')
                # take everything in the cluster assignment string except for the final token,
                # which is the object's id within the cluster.
                cl_map.setdefault(tuple(['orig'] + hierarchy[:-1]), []).append(fields[-1])

            # rename clusters and order descending in size
            desc_key = sorted(cl_map, key=lambda x: len(cl_map[x]), reverse=True)
            for n, k in enumerate(desc_key):
                cl_map[n] = np.array(cl_map.pop(k), dtype=np.int)

        return cl_map

    def _write_edges(g, parent_dir, base_name, sep=' '):
        """
        Prepare an edge-list file from the specified graph. This will be written to the
        specified parent directory, using basename.edges
        :param g: the graph
        :param parent_dir: parent directory of file
        :param base_name: file base name
        :param sep: separator within file
        :return: file name
        """
        edge_file = os.path.join(parent_dir, '{}.edges'.format(base_name))
        nx.write_edgelist(g, edge_file, data=['weight'], delimiter=sep)
        return edge_file

    assert os.path.exists(work_dir), 'supplied output path [{}] does not exist'.format(work_dir)

    base_name = 'cm_graph'
    g = to_graph(contact_map, min_len=min_len, min_sig=min_sig, norm=True, bisto=True, scale=True)

    method = method.lower()
    logger.info('Clustering contact graph using method: {}'.format(method))

    if method == 'louvain':
        cl_to_ids = louvain.cluster(g, no_iso=False, ragbag=False)
    elif method == 'mcl':
        with open(os.path.join(work_dir, 'mcl.log'), 'w+') as stdout:
            ofile = os.path.join(work_dir, '{}.mcl'.format(base_name))
            edge_file = _write_edges(g, work_dir, base_name)
            nx.write_edgelist(g, edge_file, data=['weight'])
            subprocess.check_call([app_path('external', 'mcl'),  edge_file, '--abc', '-I', '1.2', '-o', ofile],
                                  stdout=stdout, stderr=subprocess.STDOUT)
            cl_to_ids = _read_mcl(ofile)
    elif method == 'simap':
        with open(os.path.join(work_dir, 'simap.log'), 'w+') as stdout:
            ofile = os.path.join(work_dir, '{}.simap'.format(base_name))
            edge_file = _write_edges(g, work_dir, base_name)
            subprocess.check_call(['java', '-jar', app_path('external', 'simap-1.0.0.jar'), 'mdl', '-s', str(seed),
                                   '-i', '1e-5', '1e-3', '-a', '1e-5', '-g', edge_file, '-o', ofile],
                                  stdout=stdout, stderr=subprocess.STDOUT)
            cl_to_ids = _read_table(ofile)
    elif method == 'infomap':
        with open(os.path.join(work_dir, 'infomap.log'), 'w+') as stdout:
            edge_file = _write_edges(g, work_dir, base_name)
            subprocess.check_call([app_path('external', 'Infomap'), '-u', '-v', '-z', '-i', 'link-list',
                                   '-s', str(seed), '-N', '10', edge_file, work_dir],
                                  stdout=stdout, stderr=subprocess.STDOUT)
            cl_to_ids = _read_tree(os.path.join(work_dir, '{}.tree'.format(base_name)))
    elif method == 'slm':
        with open(os.path.join(work_dir, 'slm.log'), 'w+') as stdout:
            mod_func = '1'
            resolution = '2.0'
            opti_algo = '3'
            n_starts = '10'
            n_iters = '10'
            ofile = os.path.join(work_dir, '{}.slm'.format(base_name))
            verb = '1'
            edge_file = _write_edges(g, work_dir, base_name, sep='\t')
            subprocess.check_call(['java', '-jar', app_path('external', 'ModularityOptimizer.jar'), edge_file,
                                   ofile, mod_func, resolution, opti_algo, n_starts, n_iters, str(seed), verb],
                                  stdout=stdout, stderr=subprocess.STDOUT)
            cl_to_ids = _read_table(ofile, seq_col=None, cl_col=0)
    else:
        raise RuntimeError('unimplemented method: {}'.format(method))

    logger.info('Clustering using {} resulted in {} clusters'.format(method, len(cl_to_ids)))

    # standardise the results, where sequences in each cluster
    # are listed in ascending order
    clustering = {}
    for cl_id, _seqs in cl_to_ids.iteritems():
        _ord = SeqOrder.asindex(np.sort(_seqs))
        # IMPORTANT!! sequences are remapped to their gapless indices
        _seqs = contact_map.order.remap_gapless(_ord)['index']

        clustering[cl_id] = {
            'seq_ids': _seqs,
            'extent': contact_map.order.lengths()[_seqs].sum()
            # TODO add other details for clusters here
            # - residual modularity, permanence
        }

    # reestablish clusters in descending order of extent
    sorted_keys = sorted(clustering, key=lambda k: clustering[k]['extent'], reverse=True)
    clustering = {n: clustering[k] for n, k in enumerate(sorted_keys)}

    add_cluster_names(clustering)

    return clustering


def cluster_report(contact_map, clustering, source_fasta=None, is_spades=True):
    """
    For each cluster, analyze the member sequences and build a report.
    Update the clustering dictionary with this result by adding a "report" for each.

    :param contact_map: an instance of ContactMap to cluster
    :param clustering: clustering solution dictionary
    :param source_fasta: source assembly fasta (other than defined at instantiation)
    :param is_spades: if SPAdes output, we can extract coverage information from the sequence names
    """

    logger.info('Analyzing the contents of each cluster')

    seq_info = contact_map.seq_info

    if source_fasta is None:
        source_fasta = contact_map.seq_file

    # set up indexed access to the input fasta
    logger.info('Building random access index for input FASTA sequences')
    with contextlib.closing(IndexedFasta(source_fasta)) as seq_db:
        # iterate over the cluster set, in the existing order
        for cl_id, cl_info in tqdm.tqdm(clustering.iteritems(), total=len(clustering),
                                        desc='inspecting clusters'):
            _len = []
            _cov = []
            _gc = []
            for n, _seq_id in enumerate(np.sort(cl_info['seq_ids']), 1):
                # get the sequence's external name and length
                _name = seq_info[_seq_id].name
                _len.append(seq_info[_seq_id].length)
                # fetch the SeqRecord object from the input fasta
                _seq = seq_db[_name]
                _gc.append(SeqUtils.GC(_seq.seq))
                if is_spades:
                    _cov.append(float(_name.split('_')[-1]))

            if is_spades:
                report = np.array(zip(_len, _gc, _cov),
                                  dtype=[('length', np.int),
                                         ('gc', np.float),
                                         ('cov', np.float)])
            else:
                report = np.array(zip(_len, _gc),
                                  dtype=[('length', np.int),
                                         ('gc', np.float)])
            clustering[cl_id]['report'] = report


def to_graph(contact_map, norm=True, bisto=False, scale=False, extern_ids=False, min_len=None, min_sig=None):
    """
    Convert the seq_map to a undirected Networkx Graph.

    The contact map is effectively an adjacency matrix, where sequences
    are the nodes and edges weighted by the observed counts. Self-loops
    are not included by default and weights are affected by normalisation
    choices.

    :param contact_map: an instance of ContactMap to cluster
    :param norm: normalize weights by length
    :param bisto: normalise using bistochasticity
    :param scale: scale weights (max_w = 1)
    :param extern_ids: use the original external sequence identifiers for node ids
    :param min_len: override minimum sequence length, otherwise use instance's setting)
    :param min_sig: override minimum off-diagonal signal (in raw counts), otherwise use instance's setting)
    :return: graph of contigs
    """
    if extern_ids:
        _nn = lambda x: contact_map.seq_info[x].name
    else:
        _nn = lambda x: x

    if not min_len and not min_sig:
        # use the (potentially) existing mask if default criteria
        contact_map.set_primary_acceptance_mask()
    else:
        # update the acceptance mask if user has specified new criteria
        contact_map.set_primary_acceptance_mask(min_len, min_sig, update=True)

    if contact_map.processed_map is None:
        contact_map.prepare_seq_map(norm=norm, bisto=bisto)
    _map = contact_map.get_subspace(marginalise=True, flatten=False)

    logger.info('Graph will have {} nodes'.format(contact_map.order.count_accepted()))

    if not sp.isspmatrix_coo(_map):
        _map = _map.tocoo()
    scl = 1.0/_map.max() if scale else 1

    logger.debug('Building graph from edges')
    g = nx.Graph(name='contact_graph')
    for u, v, w in tqdm.tqdm(itertools.izip(_map.row, _map.col, _map.data), desc='adding edges', total=_map.nnz):
        g.add_edge(_nn(u), _nn(v), weight=w * scl)

    logger.info('Finished: {}'.format(nx.info(g).replace('\n', ' ')))

    return g


def enable_clusters(contact_map, clustering, cl_list=None, ordered_only=True, min_extent=None):
    """
    Given a clustering and list of cluster ids (or none), enable (unmask) the related sequences in
    the contact map. If a requested cluster has not been ordered, it will be dropped.

    :param contact_map: an instance of ContactMap to cluster
    :param clustering: a clustering solution to the contact map
    :param cl_list: a list of cluster ids to enable or None (all ordered clusters)
    :param ordered_only: include only clusters which have been ordered
    :param min_extent: include only clusters whose total extent is greater
    :return: the filtered list of cluster ids in ascending numerical order
    """

    # start with all clusters if unspecified
    if cl_list is None:
        cl_list = clustering.keys()

    # use instance criterion if not explicitly set
    if min_extent is None:
        min_extent = contact_map.min_extent

    if min_extent:
        cl_list = [k for k in cl_list if clustering[k]['extent'] >= min_extent]
        logger.info('Clusters passing minimum extent criterion: {}'.format(len(cl_list)))
        if len(cl_list) == 0:
            raise NoRemainingClustersException(
                'No clusters passed min_extent criterion of >= {}'.format(min_extent))

    if ordered_only:
        # drop any clusters that have not been ordered
        cl_list = [k for k in cl_list if 'order' in clustering[k]]
        logger.info('Clusters passing ordered-only criterion: {}'.format(len(cl_list)))
        if len(cl_list) == 0:
            raise NoRemainingClustersException(
                'No clusters passed ordered-only criterion')

    # impose a consistent order
    cl_list = sorted(cl_list)

    if ordered_only:
        # use the determined order and orientation
        cmb_ord = np.hstack([clustering[k]['order'] for k in cl_list])
    else:
        # arbitrary order and orientation
        cmb_ord = np.hstack([SeqOrder.asindex(clustering[k]['seq_ids']) for k in cl_list])

    if len(cmb_ord) == 0:
        raise NoRemainingClustersException('No requested cluster contained ordering information')

    logger.info('Total number of sequences in the clustering: {}'.format(len(cmb_ord)))

    # prepare the mask
    _mask = np.zeros_like(contact_map.order.mask_vector(), dtype=np.bool)
    _mask[cmb_ord['index']] = True
    _mask &= contact_map.get_primary_acceptance_mask()
    logger.info('After joining with active sequence mask map: {}'.format(_mask.sum()))
    contact_map.order.set_mask_only(_mask)
    contact_map.order.set_order_and_orientation(cmb_ord, implicit_excl=True)

    return cl_list


def plot_clusters(contact_map, fname, clustering, cl_list=None, simple=True, permute=False, max_image_size=None,
                  ordered_only=False, min_extent=None, use_taxo=False, flatten=False, **kwargs):
    """
    Plot the contact map, annotating the map with cluster names and boundaries.

    For large contact maps, block reduction can be employed to reduce the size for plotting purposes. Using
    block_reduction=2 will reduce the map dimensions by a factor of 2. Must be integer.

    :param contact_map: an instance of ContactMap to cluster
    :param fname: output file name
    :param clustering: the cluster solution
    :param cl_list: the list of cluster ids to include in plot. If none, include all ordered clusters
    :param simple: True plot seq map, False plot the extent map
    :param permute: permute the map with the present order
    :param max_image_size:  maximum allowable image size before rescale occurs
    :param ordered_only: include only clusters which have been ordered
    :param min_extent: include only clusters whose total extent is greater
    :param use_taxo: use taxonomic information within clustering, assuming it exists
    :param flatten: for tip-based, flatten matrix rather than marginalise
    :param kwargs: additional options passed to plot()
    """

    if cl_list is None:
        logger.info('Plotting heatmap of complete solution')
    else:
        logger.info('Plotting heatmap for {} specified clusters'.format(len(cl_list)))

    if simple or contact_map.bin_size is None:
        # prepare the map early as we wish to override the mask
        # which happens to be initialized in this method call
        if contact_map.processed_map is None:
            contact_map.prepare_seq_map(norm=True, bisto=True)

    # now build the list of relevant clusters and setup the associated mask
    cl_list = enable_clusters(contact_map, clustering, cl_list=cl_list, ordered_only=ordered_only,
                              min_extent=min_extent)

    if simple or contact_map.bin_size is None:
        # tick spacing simple the number of sequences in the cluster
        tick_locs = np.cumsum([0] + [len(clustering[k]['seq_ids']) for k in cl_list])
        if contact_map.is_tipbased() and flatten:
            tick_locs *= 2
    else:
        # tick spacing depends on cumulative bins for sequences in cluster
        # cumulative bin count, excluding masked sequences
        csbins = [0]
        for k in cl_list:
            # get the order records for the sequences in cluster k
            _oi = contact_map.order.order[clustering[k]['seq_ids']]
            # count the cumulative bins at each cluster for those sequences which are not masked
            csbins.append(contact_map.grouping.bins[clustering[k]['seq_ids'][_oi['mask']]].sum() + csbins[-1])
        tick_locs = np.array(csbins, dtype=np.int)

    if use_taxo:
        _labels = [clustering[cl_id]['taxon'] for cl_id in cl_list]
    else:
        _labels = [clustering[cl_id]['name'] for cl_id in cl_list]

    contact_map.plot(fname, permute=permute, simple=simple, tick_locs=tick_locs, tick_labs=_labels,
                     max_image_size=max_image_size, flatten=flatten, **kwargs)


def write_report(fname, clustering):
    """
    Create a tabular report of each cluster from a clustering report. Write the table to CSV.

    :param fname: the CSV output file name
    :param clustering: the input clustering, which contains a report
    """

    def _expect(w, x):
        """
        Weighted expectation of x with weights w. Weights do not need to be
        normalised

        :param w: weights
        :param x: variable
        :return: expectation value of x
        """
        wsum = float(w.sum())
        return np.sum(w * x) / wsum

    def _n50(x):
        """
        Calculate N50 for the given list of sequence lengths.
        :param x: a list of sequence lengths
        :return: the N50 value
        """
        x = np.sort(x)[::-1]
        return x[x.cumsum() > x.sum() / 2][0]

    df = []
    has_cov = False
    for k, v in clustering.iteritems():
        try:
            sr = v['report']

            _cl_info = [k,
                        v['name'],
                        len(v['seq_ids']),
                        v['extent'],
                        _n50(sr['length']),
                        _expect(sr['length'],
                        sr['gc']),
                        sr['gc'].mean(),
                        np.median(sr['gc']),
                        sr['gc'].std()]

            # if coverage information exists, add statistics to the table
            if 'cov' in sr.dtype.names:
                has_cov = True
                _cl_info.extend([_expect(sr['length'],
                                sr['cov']),
                                sr['cov'].mean(),
                                np.median(sr['cov']),
                                sr['cov'].std()])

            df.append(_cl_info)

        except KeyError:
            raise NoReportException(k)

    _cols = ['id', 'name', 'size', 'extent', 'n50', 'gc_expect', 'gc_mean', 'gc_median', 'gc_std']
    if has_cov:
        _cols.extend(['cov_expect', 'cov_mean', 'cov_median', 'cov_std'])

    df = pandas.DataFrame(df, columns=_cols)
    df.set_index('id', inplace=True)
    df.to_csv(fname, sep=',')


def write_mcl(contact_map, fname, clustering):
    """
    Write out the clustering solution in the format used by MCL. Each line represents a cluster
    with all members on the line show as a space-delimited list.

    :param contact_map: an instance of ContactMap to cluster
    :param fname: output file name
    :param clustering: our clustering solution
    """
    with open(fname, 'w') as outh:
        seq_info = contact_map.seq_info
        # track those sequences that were rejected during filtering
        lost = np.ones(contact_map.total_seq, dtype=np.bool)
        cl_soln = {}
        for k, v in clustering.iteritems():
            # sequence wasn't lost to filtering
            lost[v['seq_ids']] = False
            cl_soln[k] = [seq_info[ix].name for ix in np.sort(v['seq_ids'])]

        # create singleton clusters for all the lost sequences.
        # if these are left out, scoring measures aren't happy
        for n, ix in enumerate(np.argwhere(lost), len(cl_soln)):
            cl_soln[n] = [seq_info[ix[0]].name]

        clid_ascending = sorted(cl_soln.keys())
        for k in clid_ascending:
            outh.write(' '.join(cl_soln[k]))
            outh.write('\n')


def write_fasta(contact_map, output_dir, clustering, source_fasta=None, clobber=False, only_large=False):
    """
    Write out multi-fasta for all determined clusters in clustering.

    For each cluster, sequence order and orientation is as follows.
    1. for unordered clusters, sequences will be in descending nucleotide length and
       in original input orientation.
    2. for ordered clusters, sequences will appear in the prescribed order and
       orientation.

    :param contact_map: an instance of ContactMap to cluster
    :param output_dir: parent output path
    :param clustering: the clustering result, possibly also ordered
    :param source_fasta: specify a source fasta file, otherwise assume the same path as was used in parsing
    :param clobber: True overwrite files in the output path. Does not remove directories
    :param only_large: Limit output to only clusters whose extent exceedds min_extent setting
    """

    make_dir(output_dir, exist_ok=True)

    logger.info('Writing output to the path: {}'.format(output_dir))

    seq_info = contact_map.seq_info

    parent_dir = os.path.join(output_dir, 'fasta')
    make_dir(parent_dir)

    if source_fasta is None:
        source_fasta = contact_map.seq_file

    # set up indexed access to the input fasta
    with contextlib.closing(IndexedFasta(source_fasta)) as seq_db:

        # iterate over the cluster set, in the existing order
        for cl_id, cl_info in clustering.iteritems():

            if only_large and cl_info['extent'] < contact_map.min_extent:
                continue

            # Each cluster produces a multi-fasta. Sequences are not joined
            cl_path = os.path.join(parent_dir, '{}.fna'.format(cl_info['name']))

            if not clobber and os.path.exists(cl_path):
                raise IOError('Output path exists [{}] and overwriting not enabled'.format(cl_path))

            # determine the number of digits required for cluster sequence names
            try:
                num_width = max(1, int(np.ceil(np.log10(len(cl_info['seq_ids'])+1))))
            except OverflowError:
                num_width = 1

            with open(cl_path, 'w') as output_h:

                logger.debug('Writing full unordered FASTA for cluster {} to {}'.format(cl_id, cl_path))

                # iterate simply over sequence ids, while imposing ascending numerical order
                for n, _seq_id in enumerate(np.sort(cl_info['seq_ids']), 1):

                    # get the sequence's external name and length
                    _name = seq_info[_seq_id].name
                    _length = seq_info[_seq_id].length
                    # fetch the SeqRecord object from the input fasta
                    _seq = seq_db[_name]
                    # orientations are listed as unknown
                    _ori_symb = 'UNKNOWN'

                    # add a new name and description
                    _seq.id = '{0}_{1:0{2}d}'.format(cl_info['name'], n, num_width)
                    _seq.name = _seq.id
                    _seq.description = 'contig:{} ori:{} length:{}'.format(_name, _ori_symb, _length)
                    SeqIO.write(_seq, output_h, 'fasta')

            # write a separate ordered fasta as this is often a subset of all sequences
            if 'order' in cl_info:

                # Each cluster produces a multi-fasta. Sequences are not joined
                cl_path = os.path.join(parent_dir, '{}.ordered.fna'.format(cl_info['name']))

                if not clobber and os.path.exists(cl_path):
                    raise IOError('Output path exists [{}] and overwriting not enabled'.format(cl_path))

                with open(cl_path, 'w') as output_h:

                    logger.debug('Writing ordered FASTA for cluster {} to {}'.format(cl_id, cl_path))

                    # iterate over cluster members, in the determined order
                    for n, _oi in enumerate(cl_info['order'], 1):

                        # get the sequence's external name and length
                        _name = seq_info[_oi['index']].name
                        _length = seq_info[_oi['index']].length
                        # fetch the SeqRecord object from the input fasta
                        _seq = seq_db[_name]

                        # reverse complement as needed
                        if _oi['ori'] == SeqOrder.REVERSE:
                            _seq = _seq.reverse_complement()
                            _ori_symb = '-'
                        elif _oi['ori'] == SeqOrder.FORWARD:
                            _ori_symb = '+'
                        else:
                            raise UnknownOrientationStateException(_oi['ori'])

                        # add a new name and description
                        _seq.id = '{0}_{1:0{2}d}'.format(cl_info['name'], n, num_width)
                        _seq.name = _seq.id
                        _seq.description = 'contig:{} ori:{} length:{}'.format(_name, _ori_symb, _length)
                        SeqIO.write(_seq, output_h, 'fasta')
