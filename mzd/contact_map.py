import matplotlib
matplotlib.use('Agg')

from mzd import io_utils, sparse_utils
from mzd.seq_utils import *
from mzd.utils import *
from collections import OrderedDict, namedtuple
from functools import partial
from numba import jit, int64, float64, void
import Bio.SeqIO as SeqIO
import logging
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import pysam
import seaborn
import scipy.sparse as sp
import tqdm

# package logger
logger = logging.getLogger(__name__)


SeqInfo = namedtuple('SeqInfo', ['offset', 'refid', 'name', 'length', 'sites'])


def mean_selector(name):
    """
    Basic Mean functions
    """
    def geometric_mean(x, y):
        return (x*y)**0.5

    def harmonic_mean(x, y):
        return 2*x*y/(x+y)

    def arithmetic_mean(x, y):
        return 0.5*(x+y)

    try:
        mean_switcher = {
            'geometric': geometric_mean,
            'harmonic': harmonic_mean,
            'arithmetic': arithmetic_mean
        }
        return mean_switcher[name]
    except KeyError:
        raise RuntimeError('unsupported mean type [{}]'.format(name))


@jit(int64(int64[:, :], int64), nopython=True)
def find_nearest_jit(group_sites, x):
    """
    Find the nearest site from a given position on a contig.

    :param group_sites:
    :param x: query position
    :return: tuple of site and group number
    """
    ix = np.searchsorted(group_sites[:, 0], x)
    if ix == len(group_sites):
        # raise RuntimeError('find_nearest: {} didnt fit in {}'.format(x, group_sites))
        return group_sites[-1, 1]
    return group_sites[ix, 1]


@jit(void(float64[:, :, :, :], float64[:, :], float64[:], float64))
def fast_norm_tipbased_bylength(coords, data, tip_lengths, tip_size):
    """
    In-place normalisation of the sparse 4D matrix used in tip-based maps.

    As tip-based normalisation is slow for large matrices, the inner-loop has been
    moved to a Numba method.

    :param coords: the COO matrix coordinate member variable (4xN array)
    :param data:  the COO matrix data member variable (1xN array)
    :param tip_lengths: per-element min(sequence_length, tip_size)
    :param tip_size: tip size used in map
    """
    for ii in xrange(coords.shape[1]):
        i, j = coords[:2, ii]
        data[ii] *= tip_size**2 / (tip_lengths[i] * tip_lengths[j])


@jit(void(float64[:, :, :, :], float64[:, :], float64[:, :]))
def fast_norm_tipbased_bysite(coords, data, sites):
    """
    In-place normalisation of the sparse 4D matrix used in tip-based maps.

    As tip-based normalisation is slow for large matrices, the inner-loop has been
    moved to a Numba method.

    :param coords: the COO matrix coordinate member variable (4xN array)
    :param data:  the COO matrix data member variable (1xN array)
    :param sites: per-element min(sequence_length, tip_size)
    """
    for n in xrange(coords.shape[1]):
        i, j, k, l = coords[:, n]
        data[n] *= 1.0/(sites[i, k] * sites[j, l])


@jit(void(float64[:], float64[:], float64[:], float64[:]))
def fast_norm_fullseq_bysite(rows, cols, data, sites):
    """
    In-place normalisation of the scipy.coo_matrix for full sequences

    :param rows: the COO matrix coordinate member variable (4xN array)
    :param cols: the COO matrix coordinate member variable (4xN array)
    :param data:  the COO matrix data member variable (1xN array)
    :param sites: per-element min(sequence_length, tip_size)
    """
    for n in xrange(data.shape[0]):
        i = rows[n]
        j = cols[n]
        data[n] *= 1.0/(sites[i] * sites[j])


class ExtentGrouping:

    def __init__(self, seq_info, bin_size):
        self.bins = []
        self.bin_size = bin_size
        self.map = []
        self.borders = []
        self.centers = []
        self.total_bins = 0

        for n, seq in tqdm.tqdm(enumerate(seq_info), total=len(seq_info), desc='Making bins'):

            if seq.length == 0:
                raise ZeroLengthException(seq.id)

            # integer bin estimation
            num_bins = seq.length / bin_size
            if num_bins == 0:
                num_bins += 1
            # handle non-integer discrepancy by contracting/expanding all bins equally
            # the threshold between contract/expand being half a bin size
            if seq.length % bin_size != 0 and seq.length/float(bin_size) - num_bins >= 0.5:
                num_bins += 1

            edges = np.linspace(0, seq.length, num_bins+1, endpoint=True, dtype=np.int)

            self.bins.append(num_bins)

            # Per reference coordinate pairs (bin_edge, map_index)
            first_bin = self.total_bins
            last_bin = first_bin + num_bins
            self.map.append(np.vstack((edges[1:], np.arange(first_bin, last_bin))).T)
            self.borders.append(np.array([first_bin, last_bin], dtype=np.int))

            self.total_bins += num_bins

            c_nk = edges[:-1] + 0.5*(edges[1] - edges[0]) - 0.5*seq.length
            self.centers.append(c_nk.reshape((1, len(c_nk))))
            # logger.debug('{}: {} bins'.format(n, num_bins))

        self.bins = np.array(self.bins)


class SeqOrder:

    FORWARD = 1
    REVERSE = -1

    ACCEPTED = True
    EXCLUDED = False

    STRUCT_TYPE = np.dtype([('pos', np.int32), ('ori', np.int8), ('mask', np.bool), ('length', np.int32)])
    INDEX_TYPE = np.dtype([('index', np.int32), ('ori', np.int8)])

    def __init__(self, seq_info):
        """
        Initial order is determined by the order of supplied sequence information dictionary. Sequences
        are given surrogate ids using consecutive integers. Member functions expect surrogate ids
        not original names.

        The class also retains orientation and masking state. Orientation defines whether a sequence
        should be in its original direction (as read in) (1) or reverse complemented (-1).

        Masking state defines whether a input sequence shall been excluded from further consideration.
        (accepted=1, excluded=0)

        :param seq_info: sequence information dictionary
        """
        self._positions = None
        _ord = np.arange(len(seq_info), dtype=np.int32)
        self.order = np.array(
            [(_ord[i], SeqOrder.FORWARD, SeqOrder.ACCEPTED, seq_info[i].length) for i in xrange(len(_ord))],
            dtype=SeqOrder.STRUCT_TYPE)

        self._update_positions()

    @staticmethod
    def asindex(_ord):
        """
        Convert a simple list or ndarray of indices, to INDEX_TYPE array with default forward orientation.

        :param _ord: list/ndarray of indices
        :return: INDEX_TYPE array
        """
        assert isinstance(_ord, (list, np.ndarray)), 'input must be a list or ndarray'
        return np.array(zip(_ord, np.ones_like(_ord, dtype=np.bool)), dtype=SeqOrder.INDEX_TYPE)

    def _update_positions(self):
        """
        An optimisation, whenever the positional state changes, this method must be called to
        maintain the current state in a separate array. This avoids unnecessary recalculation
        overhead.
        """
        # Masked sequences last, then by current position.
        sorted_indices = np.lexsort([self.order['pos'], ~self.order['mask']])
        for n, i in enumerate(sorted_indices):
            self.order[i]['pos'] = n
        self._positions = np.argsort(self.order['pos'])

    def remap_gapless(self, gapless_indices):
        """
        Recover the original, potentially sparse (gapped) indices from a dense (gapless) set
        of indices. Gaps originate from sequences being masked in the order. External tools
        often expect and return dense indices. When submitting changes to the current order
        state, it is important to first apply this method and reintroduce any gaps.

        Both a list/array of indices or a INDEX_TYPE array can be passed.

        :param gapless_indices: dense list of indices or an ndarray of type INDEX_TYPE
        :return: remappped indices with gaps (of a similar type to input)
        """
        # not as yet verified but this method is being replaced by the 50x faster numpy
        # alternative below. The slowless shows for large problems and repeated calls.
        # we ~could~ go further and maintain the shift array but this will require
        # consistent and respectful (fragile) use of mutator methods and not direct access on mask
        # or an observer.

        # the accumulated shifts due to masked sequences (the gaps).
        # we remove the masked sequences to make this array gapless
        shift = np.cumsum(~self.order['mask'])[self.order['mask']]

        # now reintroduce the gaps to the gapless representation supplied

        # handle our local type
        if isinstance(gapless_indices, np.ndarray) and gapless_indices.dtype == SeqOrder.INDEX_TYPE:
            remapped = []
            for oi in gapless_indices:
                remapped.append((oi['index'] + shift[oi['index']], oi['ori']))
            return np.array(remapped, dtype=SeqOrder.INDEX_TYPE)

        # handle plain collection
        else:
            remapped = []
            for oi in gapless_indices:
                remapped.append(oi + shift[oi])
            return np.array(remapped)

    def accepted_positions(self, copy=True):
        """
        The current positional order of only those sequences which have not been excluded by the mask.
        :param copy: return a copy
        :return: all accepted positons in order of index
        """
        return self.all_positions(copy=copy)[:self.count_accepted()]

    def all_positions(self, copy=True):
        """
        The current positional order of all sequences. Internal logic relegates masked sequences to always come
        last and ascending surrogate id order.

        :param copy: return a copy of the positions
        :return: all positions in order of index, masked or not.
        """
        if copy:
            _p = self._positions.copy()
        else:
            _p = self._positions
        return _p

    @staticmethod
    def double_order(_ord):
        """
        For doublet maps, the stored order must be re-expanded to reference the larger (2x) map.

        :param _ord:
        :return: expanded order
        """
        return np.array([[2*oi, 2*oi+1] for oi in _ord]).ravel()

    def gapless_positions(self):
        """
        A dense index range representing the current positional order without masked sequences. Therefore
        the returned array does not contain surrogate ids, but rather the relative positions of unmasked
        sequences, when all masked sequences have been discarded.

        :return: a dense index range of positional order, once all masked sequences have been discarded.
        """
        # accumulated shift from gaps
        gap_shift = np.cumsum(~self.order['mask'])
        # just unmasked sequences
        _p = np.argsort(self.order['pos'])
        _p = _p[:self.count_accepted()]
        # removing gaps leads to a dense range of indices
        _p -= gap_shift[_p]
        return _p

    def set_mask_only(self, _mask):
        """
        Set the mask state of all sequences, where indices in the mask map to
        sequence surrogate ids.

        :param _mask: mask array or list, boolean or 0/1 valued
        """
        _mask = np.asarray(_mask, dtype=np.bool)
        assert len(_mask) == len(self.order), 'supplied mask must be the same length as existing order'
        assert np.all((_mask == SeqOrder.ACCEPTED) | (_mask == SeqOrder.EXCLUDED)), \
            'new mask must be {} or {}'.format(SeqOrder.ACCEPTED, SeqOrder.EXCLUDED)

        # assign mask
        self.order['mask'] = _mask
        self._update_positions()

    def set_order_only(self, _ord, implicit_excl=False):
        """
        Convenience method to set the order using a list or 1D ndarray. Orientations will
        be assumed as all forward (+1).

        :param _ord: a list or ndarray of surrogate ids
        :param implicit_excl: implicitly extend the order to include unmentioned excluded sequences.
        """
        assert isinstance(_ord, (list, np.ndarray)), 'Wrong type supplied, order must be a list or ndarray'
        if isinstance(_ord, np.ndarray):
            _ord = np.ravel(_ord)
            assert np.ndim(_ord) == 1, 'orders as numpy arrays must be 1-dimensional'
        # augment the order to include default orientations
        _ord = SeqOrder.asindex(_ord)
        self.set_order_and_orientation(_ord, implicit_excl=implicit_excl)

    def set_order_and_orientation(self, _ord, implicit_excl=False):
        """
        Set only the order, while ignoring orientation. An ordering is defined
        as a 1D array of the structured type INDEX_TYPE, where elements are the
        position and orientation of each indices.

        NOTE: This definition can be the opposite of what is returned by some
        ordering methods, and np.argsort(_v) should inverse the relation.

        NOTE: If the order includes only active sequences, setting implicit_excl=True
        the method will implicitly assume unmentioned ids are those currently
        masked. An exception is raised if a masked sequence is included in the order.

        :param _ord: 1d ordering
        :param implicit_excl: implicitly extend the order to include unmentioned excluded sequences.
        """
        assert _ord.dtype == SeqOrder.INDEX_TYPE, 'Wrong type supplied, _ord should be of INDEX_TYPE'

        if len(_ord) < len(self.order):
            # some sanity checks
            assert implicit_excl, 'Use implicit_excl=True for automatic handling ' \
                                  'of orders only mentioning accepted sequences'
            assert len(_ord) == self.count_accepted(), 'new order must mention all ' \
                                                       'currently accepted sequences'
            # those surrogate ids mentioned in the order
            mentioned = set(_ord['index'])
            assert len(mentioned & set(self.excluded())) == 0, 'new order and excluded must not ' \
                                                               'overlap when using implicit assignment'
            assert len(mentioned ^ set(self.accepted())) == 0, 'incomplete new order supplied,' \
                                                               'missing accepted ids'
            # assign the new orders
            self.order['pos'][_ord['index']] = np.arange(len(_ord), dtype=np.int32)
            self.order['ori'][_ord['index']] = _ord['ori']
            # mask remaining, unmentioned indices
            _mask = np.zeros_like(self.mask_vector(), dtype=np.bool)
            _mask[_ord['index']] = True
            self.set_mask_only(_mask)
        else:
            # just a simple complete order update
            assert len(_ord) == len(self.order), 'new order was a different length'
            assert len(set(_ord['index']) ^ set(self.accepted())) == 0, 'incomplete new order supplied,' \
                                                                        'missing accepted ids'
            self.order['pos'][_ord['index']] = np.arange(len(_ord), dtype=np.int32)
            self.order['ori'][_ord['index']] = _ord['ori']

        self._update_positions()

    def accepted_order(self):
        """
        :return: an INDEX_TYPE array of the order and orientation of the currently accepted sequences.
        """
        idx = np.where(self.order['mask'])
        ori = np.ones(self.count_accepted(), dtype=np.int)
        return np.array(zip(idx, ori), dtype=SeqOrder.INDEX_TYPE)

    def mask_vector(self):
        """
        :return: the current mask vector
        """
        return self.order['mask']

    def mask(self, _id):
        """
        Mask an individual sequence by its surrogate id

        :param _id: the surrogate id of sequence
        """
        self.order[_id]['mask'] = False
        self._update_positions()

    def count_accepted(self):
        """
        :return: the current number of accepted (unmasked) sequences
        """
        return self.order['mask'].sum()

    def count_excluded(self):
        """
        :return: the current number of excluded (masked) sequences
        """
        return len(self.order) - self.count_accepted()

    def accepted(self):
        """
        :return: the list surrogate ids for currently accepted sequences
        """
        return np.where(self.order['mask'])[0]

    def excluded(self):
        """
        :return: the list surrogate ids for currently excluded sequences
        """
        return np.where(~self.order['mask'])[0]

    def flip(self, _id):
        """
        Flip the orientation of the sequence

        :param _id: the surrogate id of sequence
        """
        self.order[_id]['ori'] *= -1

    def lengths(self, exclude_masked=False):
        # type: (bool) -> np.ndarray
        """
        Sequence lengths

        :param exclude_masked: True include only umasked sequencces
        :return: the lengths of sequences
        """
        if exclude_masked:
            return self.order['length'][self.order['mask']]
        else:
            return self.order['length']

    def shuffle(self):
        """
        Randomize order
        """
        np.random.shuffle(self.order['pos'])
        self._update_positions()

    def before(self, a, b):
        """
        Test if a comes before another sequence in order.

        :param a: surrogate id of sequence a
        :param b: surrogate id of sequence b
        :return: True if a comes before b
        """
        assert a != b, 'Surrogate ids must be different'
        return self.order['pos'][a] < self.order['pos'][b]

    def intervening(self, a, b):
        """
        For the current order, calculate the length of intervening
        sequences between sequence a and sequence b.

        :param a: surrogate id of sequence a
        :param b: surrogate id of sequence b
        :return: total length of sequences between a and b.
        """
        assert a != b, 'Surrogate ids must be different'

        pa = self.order['pos'][a]
        pb = self.order['pos'][b]
        if pa > pb:
            pa, pb = pb, pa
        inter_ix = self._positions[pa+1:pb]
        return np.sum(self.order['length'][inter_ix])


class ContactMap:

    def __init__(self, bam_file, enzymes, seq_file, min_insert, min_mapq=0, min_len=0, min_sig=1, min_extent=0,
                 min_size=0, max_fold=None, random_seed=None, strong=None, bin_size=None, tip_size=None,
                 precount=False):

        self.strong = strong
        self.bam_file = bam_file
        self.bin_size = bin_size
        self.min_mapq = min_mapq
        self.min_insert = min_insert
        self.min_len = min_len
        self.min_sig = min_sig
        self.min_extent = min_extent
        self.min_size = min_size
        self.max_fold = max_fold
        self.random_state = np.random.RandomState(random_seed)
        self.strong = strong
        self.seq_info = []
        self.seq_map = None
        self.seq_file = seq_file
        self.grouping = None
        self.extent_map = None
        self.order = None
        self.tip_size = tip_size
        self.precount = precount
        self.total_reads = None
        self.cov_info = None
        self.processed_map = None
        self.primary_acceptance_mask = None
        self.bisto_scale = None
        self.seq_analyzer = None
        self.enzymes = enzymes

        # build a dictionary of sites in each reference first
        fasta_info = {}
        with io_utils.open_input(seq_file) as multi_fasta:
            # prepare the site counter for the given experimental conditions
            site_counter = SiteCounter(enzymes, tip_size, is_linear=True)
            # get an estimate of sequences for progress
            fasta_count = count_fasta_sequences(seq_file)
            for seqrec in tqdm.tqdm(SeqIO.parse(multi_fasta, 'fasta'), total=fasta_count, desc='Analyzing sites'):
                if len(seqrec) < min_len:
                    continue
                fasta_info[seqrec.id] = {'sites': site_counter.count_sites(seqrec.seq),
                                         'length': len(seqrec)}

        # now parse the header information from bam file
        with pysam.AlignmentFile(bam_file, 'rb') as bam:

            # test that BAM file is the correct sort order
            if 'SO' not in bam.header['HD'] or bam.header['HD']['SO'] != 'queryname':
                raise IOError('BAM file must be sorted by read name')

            # determine the set of active sequences
            # where the first filtration step is by length
            ref_count = {'seq_missing': 0, 'too_short': 0}
            offset = 0
            logger.info('Reading sequences...')
            for n, (rname, rlen) in enumerate(zip(bam.references, bam.lengths)):

                # minimum length threshold
                if rlen < min_len:
                    ref_count['too_short'] += 1
                    continue

                try:
                    fa = fasta_info[rname]
                except KeyError:
                    logger.info('Sequence: "{}" was not present in reference fasta'.format(rname))
                    ref_count['seq_missing'] += 1
                    continue

                assert fa['length'] == rlen, \
                    'Sequence lengths in {} do not agree: bam {} fasta {}'.format(rname, fa['length'], rlen)

                self.seq_info.append(SeqInfo(offset, n, rname, rlen, fa['sites']))

                offset += rlen

            # total extent covered
            self.total_len = offset
            self.total_seq = len(self.seq_info)

            # all sequences begin as active in mask
            # when filtered, mask[i] = False
            self.current_mask = np.ones(self.total_seq, dtype=np.bool)

            if self.total_seq == 0:
                logger.info('No sequences in BAM found in FASTA')
                raise ParsingError('No sequences in BAM found in FASTA')

            logger.info('Accepted {} sequences covering {} bp'.format(self.total_seq, self.total_len))
            logger.info('References excluded: {}'.format(ref_count))

            if self.bin_size:
                logger.info('Determining bins...')
                self.grouping = ExtentGrouping(self.seq_info, self.bin_size)

            logger.info('Counting reads in bam file...')

            if self.precount:
                self.total_reads = bam.count(until_eof=True)
                logger.info('BAM file contains {0} alignments'.format(self.total_reads))
            else:
                logger.info('Skipping pre-count of BAM file, no ETA will be offered')

            # initialise the order
            self.order = SeqOrder(self.seq_info)

            # accumulate
            self._bin_map(bam)

            # create an initial acceptance mask
            self.set_primary_acceptance_mask()

    def _bin_map(self, bam):
        """
        Accumulate read-pair observations from the supplied BAM file.
        Maps are initialized here. Logical control is achieved through initialisation of the
        ContactMap instance, rather than supplying this function arguments.

        :param bam: this instance's open bam file.
        """
        import tqdm

        def _simple_match(r):
            return not r.is_unmapped and r.mapq >= _mapq

        def _strong_match(r):
            cig = r.cigartuples[-1] if r.is_reverse else r.cigartuples[0]
            return (r.mapping_quality >= _mapq and
                    not r.is_secondary and
                    not r.is_supplementary and
                    cig[0] == 0 and cig[1] >= self.strong)

        # set-up match call
        _matcher = _strong_match if self.strong else _simple_match

        def _on_tip_withlocs(p1, p2, l1, l2, _tip_size):
            tailhead_mat = np.zeros((2, 2), dtype=np.uint32)
            i = None
            j = None

            # contig1 tips won't overlap
            if l1 > 2 * _tip_size:
                if p1 < _tip_size:
                    i = 0
                elif p1 > l1 - _tip_size:
                    i = 1

            # contig1 tips will overlap
            else:
                # assign to whichever end is closest
                if p1 < l1 - p1:
                    i = 0
                elif l1 - p1 < p1:
                    i = 1

            # only bother with second tip assignment if the first was ok
            if i is not None:

                # contig2 tips won't overlap
                if l2 > 2 * _tip_size:
                    if p2 < _tip_size:
                        j = 0
                    elif p2 > l2 - _tip_size:
                        j = 1

                # contig2 tips will overlap
                else:
                    # assign to whichever end is closest
                    if p2 < l2 - p2:
                        j = 0
                    elif l2 - p2 < p2:
                        j = 1

            tailhead_mat[i, j] = 1
            return i is not None and j is not None, tailhead_mat

        def _always_true(*args):
            return True, 1

        _on_tip = _always_true if not self.is_tipbased() else _on_tip_withlocs

        # initialise a sparse matrix for accumulating the map
        if not self.is_tipbased():
            # just a basic NxN sparse array for normal whole-sequence binning
            _seq_map = sparse_utils.Sparse2DAccumulator(self.total_seq)
        else:
            # each tip is tracked separately resulting in the single count becoming a 2x2 interaction matrix.
            # therefore the tensor has dimension NxNx2x2
            _seq_map = sparse_utils.Sparse4DAccumulator(self.total_seq)

        # if binning also requested, initialise another sparse matrix
        if self.bin_size:
            logger.info('Initialising contact map of {0}x{0} fragment bins, '
                        'representing {1} bp over {2} sequences'.format(self.grouping.total_bins,
                                                                        self.total_len, self.total_seq))
            _extent_map = sparse_utils.Sparse2DAccumulator(self.grouping.total_bins)
            _grouping_map = self.grouping.map
        else:
            _grouping_map = None
            _extent_map = None

        with tqdm.tqdm(total=self.total_reads) as pbar:

            # locals for read filtering
            _min_sep = self.min_insert
            _mapq = self.min_mapq

            _idx = self.make_reverse_index('refid')

            # locals for tip checking
            _len = bam.lengths
            _tip_size = self.tip_size

            counts = OrderedDict({
                'accepted': 0,
                'not_tip': 0,
                'short_insert': 0,
                'ref_excluded': 0,
                'median_excluded': 0,
                'end_buffered': 0,
                'poor_match': 0})

            bam.reset()
            bam_iter = bam.fetch(until_eof=True)
            while True:

                try:
                    r1 = bam_iter.next()
                    pbar.update()
                    while True:
                        # read records until we get a pair
                        r2 = bam_iter.next()
                        pbar.update()
                        if r1.query_name == r2.query_name:
                            break
                        r1 = r2
                except StopIteration:
                    break

                if r1.reference_id not in _idx or r2.reference_id not in _idx:
                    counts['ref_excluded'] += 1
                    continue

                if not _matcher(r1) or not _matcher(r2):
                    counts['poor_match'] += 1
                    continue

                # if no_go[r1.reference_id][r1.reference_start:r1.reference_end] or \
                #         no_go[r2.reference_id][r2.reference_start:r2.reference_end]:
                #     counts['median_excluded'] += 1
                #     continue

                if r1.is_read2:
                    r1, r2 = r2, r1

                # # accept only inward facing read pairs
                # if not ((r1.is_reverse and not r2.is_reverse) or (not r1.is_reverse and r2.is_reverse)):
                #     # counts['skipped'] += 1
                #     continue

                # use 5-prime base depending on orientation
                r1pos = r1.pos if not r1.is_reverse else r1.pos + r1.alen
                r2pos = r2.pos if not r2.is_reverse else r2.pos + r2.alen

                # filter inserts deemed "short" which tend to be heavily WGS signal
                if _min_sep and r1.is_proper_pair:
                    ins_len = r2.pos - r1.pos
                    if ins_len < _min_sep:
                        counts['short_insert'] += 1
                        continue

                # get reference lengths
                l1 = _len[r1.reference_id]
                l2 = _len[r2.reference_id]

                # get internal indices
                ix1 = _idx[r1.reference_id]
                ix2 = _idx[r2.reference_id]

                # maintain just a half-matrix
                if ix2 < ix1:
                    ix1, ix2 = ix2, ix1
                    r1pos, r2pos = r2pos, r1pos
                    l1, l2 = l2, l1

                if _extent_map:
                    b1 = find_nearest_jit(_grouping_map[ix1], r1pos)
                    b2 = find_nearest_jit(_grouping_map[ix2], r2pos)

                    # maintain half-matrix
                    if b1 > b2:
                        b1, b2 = b2, b1

                    # tally all mapped reads for binned map, not just those considered in tips
                    _extent_map[b1, b2] += 1

                # for seq-map, we may reject reads outside of a defined tip region
                tip_info = _on_tip(r1pos, r2pos, l1, l2, _tip_size)
                if not tip_info[0]:
                    counts['not_tip'] += 1
                    continue

                counts['accepted'] += 1

                _seq_map[ix1, ix2] += tip_info[1]

        # default to always making matrices symmetric
        if self.bin_size:
            self.extent_map = _extent_map.get_coo()
            del _extent_map

        self.seq_map = _seq_map.get_coo()
        del _seq_map

        logger.info('Pair accounting: {}'.format(counts))
        logger.info('Total extent map weight {}'.format(self.map_weight()))

    @staticmethod
    def get_fields():
        """
        :return: the list of fields used in seq_info dict.
        """
        return SeqInfo._fields

    def make_reverse_index(self, field_name):
        """
        Make a reverse look-up (dict) from the chosen field in seq_info to the internal index value
        of the given sequence. Non-unique fields will raise an exception.

        :param field_name: the seq_info field to use as the reverse.
        :return: internal array index of the sequence
        """
        rev_idx = {}
        for n, seq in enumerate(self.seq_info):
            fv = getattr(seq, field_name)
            if fv in rev_idx:
                raise RuntimeError('field contains non-unique entries, a 1-1 mapping cannot be made')
            rev_idx[fv] = n
        return rev_idx

    def map_weight(self):
        """
        :return: the total map weight (sum ij)
        """
        return self.seq_map.sum()

    def is_empty(self):
        """
        :return: True if the map has zero weight
        """
        return self.map_weight() == 0

    def is_tipbased(self):
        """
        :return: True if the seq_map is a tip-based 4D tensor
        """
        return self.tip_size is not None

    def get_primary_acceptance_mask(self):
        assert self.primary_acceptance_mask is not None, 'Primary acceptance mask has not be initialized'
        return self.primary_acceptance_mask.copy()

    def set_primary_acceptance_mask(self, min_len=None, min_sig=None, max_fold=None, update=False):
        """
        Determine and set the filter mask using the specified constraints across the entire
        contact map. The mask is True when a sequence is considered acceptable wrt to the
        constraints. The mask is also returned by the function for convenience.

        :param min_len: override instance value for minimum sequence length
        :param min_sig: override instance value for minimum off-diagonal signal (counts)
        :param max_fold: maximum locally-measured fold-coverage to permit
        :param update: replace the current primary mask if it exists
        :return: an acceptance mask over the entire contact map
        """
        assert max_fold is None, 'Filtering on max_fold is currently disabled'

        # If parameter based critiera were unset, use instance member values set at instantiation time
        if not min_len:
            min_len = self.min_len
        if not min_sig:
            min_sig = self.min_sig

        assert min_len, 'Filtering criteria min_len is None'
        assert min_sig, 'Filtering criteria min_sig is None'

        logger.debug('Setting primary acceptance mask with '
                     'filtering criterion min_len: {} min_sig: {}'.format(min_len, min_sig))

        # simply return the current mask if it has already been determined
        # and an update is not requested
        if not update and self.primary_acceptance_mask is not None:
            logger.debug('Using existing mask')
            return self.get_primary_acceptance_mask()

        acceptance_mask = np.ones(self.total_seq, dtype=np.bool)

        # mask for sequences shorter than limit
        _mask = self.order.lengths() >= min_len
        logger.debug('Minimum length threshold removing: {}'.format(self.total_seq - _mask.sum()))
        acceptance_mask &= _mask

        # mask for sequences weaker than limit
        if self.is_tipbased():
            signal = sparse_utils.max_offdiag_4d(self.seq_map)
        else:
            signal = sparse_utils.max_offdiag(self.seq_map)
        _mask = signal >= min_sig
        logger.debug('Minimum signal threshold removing: {}'.format(self.total_seq - _mask.sum()))
        acceptance_mask &= _mask

        # retain the union of all masks.
        self.primary_acceptance_mask = acceptance_mask

        logger.debug('Accepted sequences: {}'.format(self.primary_acceptance_mask.sum()))

        return self.get_primary_acceptance_mask()

    def prepare_seq_map(self, norm=True, bisto=False, mean_type='geometric'):
        """
        Prepare the sequence map (seq_map) by application of various filters and normalisations.

        :param norm: normalisation by sequence lengths
        :param bisto: make the output matrix bistochastic
        :param mean_type: when performing normalisation, use "geometric, harmonic or arithmetic" mean.
        """

        logger.info('Preparing sequence map with full dimensions: {}'.format(self.seq_map.shape))

        _mask = self.get_primary_acceptance_mask()

        self.order.set_mask_only(_mask)

        if self.order.count_accepted() < 1:
            raise NoneAcceptedException()

        _map = self.seq_map.astype(np.float)

        # apply length normalisation if requested
        if norm:
            _map = self._norm_seq(_map, self.is_tipbased(), mean_type=mean_type, use_sites=True)
            logger.debug('Map normalized')

        # make map bistochastic if requested
        if bisto:
            # TODO balancing may be better done after compression
            _map, scl = self._bisto_seq(_map)
            # retain the scale factors
            self.bisto_scale = scl
            logger.debug('Map balanced')

        # cache the results for optional quick access
        self.processed_map = _map

    def get_subspace(self, permute=False, external_mask=None, marginalise=False, flatten=True,
                     dtype=np.float):
        """
        Using an already normalized full seq_map, return a subspace as indicated by an external
        mask or if none is supplied, the full map without filtered elements.

        The supplied external mask must refer to all sequences in the map.

        :param permute: reorder the map with the current ordering state
        :param external_mask: an external mask to combine with the existing primary mask
        :param marginalise: Assuming 4D NxNx2x2 tensor, sum 2x2 elements to become a 2D NxN
        :param flatten: convert a NxNx2x2 tensor to a 2Nx2N matrix
        :param dtype: return map with specific element type
        :return: subspace map
        """
        assert (not marginalise and not flatten) or np.logical_xor(marginalise, flatten), \
            'marginalise and flatten are mutually exclusive'

        # starting with the normalized map
        _map = self.processed_map.astype(dtype)

        # from a union of the sequence filter and external mask
        if external_mask is not None:
            _mask = self.get_primary_acceptance_mask()
            logger.info('Beginning with sequences after primary filtering: {}'.format(_mask.sum()))
            _mask &= external_mask
            logger.info('Active sequences after applying external mask: {}'.format(_mask.sum()))
            self.order.set_mask_only(_mask)

        # remove masked sequences from the map
        if self.order.count_accepted() < self.total_seq:
            if self.is_tipbased():
                _map = sparse_utils.compress_4d(_map, self.order.mask_vector())
            else:
                _map = sparse_utils.compress(_map.tocoo(), self.order.mask_vector())
            logger.info('After removing filtered sequences map dimensions: {}'.format(_map.shape))

        # convert tip-based tensor to other forms
        if self.is_tipbased():
            if marginalise:
                logger.debug('Marginalising NxNx2x2 tensor to NxN matrix')
                # sum counts of the 2x2 confusion matrices into 1 value
                _map = _map.sum(axis=(2, 3)).to_scipy_sparse()
            elif flatten:
                logger.debug('Flattening NxNx2x2 tensor to 2Nx2N matrix')
                # convert the 4D map into a 2Nx2N 2D map.
                _map = sparse_utils.flatten_tensor_4d(_map)

        if permute:
            _map = self._reorder_seq(_map, flatten=flatten)
            logger.debug('Map reordered')

        return _map

    def get_extent_map(self, norm=True, bisto=False, permute=False, mean_type='geometric'):
        """
        Return the extent map after applying specified processing steps. Masked sequences are always removed.

        :param norm: sequence length normalisation
        :param bisto: make map bistochastic
        :param permute: permute the map using current order
        :param mean_type: length normalisation mean (geometric, harmonic, arithmetic)
        :return: processed extent map
        """

        logger.info('Preparing extent map with fill dimensions: {}'.format(self.extent_map.shape))

        _map = self.extent_map.astype(np.float)

        # apply length normalisation if requested
        if norm:
            _map = self._norm_extent(_map, mean_type)
            logger.debug('Map normalized')

        # if there are sequences to mask, remove them from the map
        if self.order.count_accepted() < self.total_seq:
            _map = self._compress_extent(_map)
            logger.info('After removing filtered sequences map dimensions: {}'.format(_map.shape))

        # make map bistochastic if requested
        if bisto:
            _map, scl = sparse_utils.kr_biostochastic(_map)
            logger.debug('Map balanced')

        # reorder using current order state
        if permute:
            _map = self._reorder_extent(_map)
            logger.debug('Map reordered')

        return _map

    def extent_to_seq(self):
        """
        Convert the extent map into a single-pixel per sequence "seq_map". This method
        is useful when only a tip based seq_map has been produced, and an analysis would be
        better done on a full accounting of mapping interactions across each sequences full
        extent.

        :return: a seq_map representing all counts across each sequence
        """
        m = self.extent_map.tocsr()
        m_out = sparse_utils.Sparse2DAccumulator(self.total_seq)
        cbins = np.cumsum(self.grouping.bins)
        a0 = 0
        for i in xrange(len(self.grouping.bins)):
            a1 = cbins[i]
            # sacrifice memory for significant speed up slicing below
            row_i = m[a0:a1, :].todense()
            b0 = 0
            for j in xrange(i, len(self.grouping.bins)):
                b1 = cbins[j]
                mij = row_i[:, b0:b1].sum()
                if mij == 0:
                    continue
                m_out[i, j] = int(mij)
                b0 = b1
            a0 = a1
        return m_out.get_coo()

    def _reorder_seq(self, _map, flatten=False):
        """
        Reorder a simple sequence map using the supplied map.

        :param _map: the map to reorder
        :param flatten: tip-based tensor converted to 2Nx2N matrix, otherwise the assumption is marginalisation
        :return: ordered map
        """
        assert sp.isspmatrix(_map), 'reordering expects a sparse matrix type'

        _order = self.order.gapless_positions()
        if self.is_tipbased() and flatten:
            _order = SeqOrder.double_order(_order)

        assert _map.shape[0] == _order.shape[0], 'supplied map and unmasked order are different sizes'
        p = sp.lil_matrix(_map.shape)
        for i in xrange(len(_order)):
            p[i, _order[i]] = 1.
        p = p.tocsr()
        return p.dot(_map.tocsr()).dot(p.T)

    def _bisto_seq(self, _map):
        """
        Make a contact map bistochastic. This is another form of normslisation. Automatically
        handles 2D and 4D maps.

        :param _map: a map to balance (make bistochastic)
        :return: the balanced map
        """
        logger.debug('Balancing contact map')

        if self.is_tipbased():
            _map, scl = sparse_utils.kr_biostochastic_4d(_map)
        else:
            _map, scl = sparse_utils.kr_biostochastic(_map)
        return _map, scl

    def _get_sites(self):
        _sites = np.array([si.sites for si in self.seq_info], dtype=np.float)
        # all sequences are assumed to have a minimum of 1 site -- even if not observed
        # TODO test whether it would be more accurate to assume that all sequences are under counted by 1.
        _sites[np.where(_sites == 0)] = 1
        return _sites

    def _norm_seq(self, _map, tip_based, use_sites=True, mean_type='geometric'):
        """
        Normalise a simple sequence map in place by the geometric mean of interacting contig pairs lengths.
        The map is assumed to be in starting order.

        :param _map: the target map to apply normalisation
        :param tip_based: treat the supplied map as a tip-based tensor
        :param use_sites: normalise matrix counts using observed sites, otherwise normalise
        using sequence lengths as a proxy
        :param mean_type: for length normalisation, choice of mean (harmonic, geometric, arithmetic)
        :return: normalized map
        """
        if use_sites:
            logger.debug('Doing site based normalisation')
            _sites = self._get_sites()
            _map = _map.astype(np.float)
            if tip_based:
                fast_norm_tipbased_bysite(_map.coords, _map.data, _sites)
            else:
                fast_norm_fullseq_bysite(_map.row, _map.col, _map.data, _sites)

        else:
            logger.debug('Doing length based normalisation')
            if tip_based:
                _tip_lengths = np.minimum(self.tip_size, self.order.lengths()).astype(np.float)
                fast_norm_tipbased_bylength(_map.coords, _map.data, _tip_lengths, self.tip_size)
            else:
                _mean_func = mean_selector(mean_type)
                _len = self.order.lengths().astype(np.float)
                _map = _map.tolil().astype(np.float)
                for i in xrange(_map.shape[0]):
                    _map[i, :] /= np.fromiter((1e-3 * _mean_func(_len[i],  _len[j])
                                               for j in xrange(_map.shape[0])), dtype=np.float)
                _map = _map.tocsr()

        return _map

    def _norm_extent(self, _map, mean_type='geometric'):
        """
        Normalise a extent map in place by the geometric mean of interacting contig pairs lengths.

        :return: a normalized extent map in lil_matrix format
        """
        assert sp.isspmatrix(_map), 'Extent matrix is not a scipy matrix type'

        if _map.dtype not in {np.float, float}:
            _map = _map.astype(np.float)
        if not sp.isspmatrix_lil(_map):
            _map = _map.tolil()

        _mean_func = mean_selector(mean_type)
        _len = self.order.lengths().astype(np.float)
        _cbins = np.cumsum(self.grouping.bins)
        for row_i, col_dat in enumerate(_map.rows):
            i = np.searchsorted(_cbins, row_i, side='right')
            wi = np.fromiter((1e-3 * _mean_func(_len[i], _len[j])
                              for j in np.searchsorted(_cbins, col_dat, side='right')), dtype=np.float)
            _map.data[row_i] /= wi
        return _map

    def _reorder_extent(self, _map):
        """
        Reorder the extent map using current order.

        :return: sparse CSR format permutation of the given map
        """
        _order = self.order.gapless_positions()
        _bins = self.grouping.bins[self.order.mask_vector()]
        _ori = self.order.order['ori'][np.argsort(self.order.order['pos'])]

        # create a permutation matrix
        p = sp.lil_matrix(_map.shape)
        _shuf_bins = _bins[_order]
        for i, oi in enumerate(_order):
            j_off = _bins[:oi].sum()
            i_off = _shuf_bins[:i].sum()
            if _ori[i] > 0:
                for k in xrange(_bins[oi]):
                    p[i_off+k, j_off+k] = 1
            else:
                # rot90 those with reverse orientation
                _nb = _bins[oi]
                for k in xrange(_nb):
                    p[i_off+_nb-(k+1), j_off+k] = 1

        # permute the extent_map
        p = p.tocsr()
        return p.dot(_map.tocsr()).dot(p.T)

    def _compress_extent(self, _map):
        """
        Compress the extent map for each sequence that is presently masked. This will eliminate
        all bins which pertain to a given masked sequence.

        :return: a scipy.sparse.coo_matrix pertaining to only the unmasked sequences.
        """
        assert sp.isspmatrix(_map), 'Extent matrix is not a scipy sparse matrix type'
        if not sp.isspmatrix_coo(_map):
            _map = _map.tocoo()

        _order = self.order.order
        _bins = self.grouping.bins

        # build a list of every accepted element.
        # TODO this could be done as below, without the memory requirements of realising all elements
        s = 0
        accept_bins = []
        # accept_index = set(np.where(_mask)[0])
        for i in xrange(len(_order)):
            # if i in accept_index:
            if _order[i]['mask']:
                accept_bins.extend([j+s for j in xrange(_bins[i])])
            s += _bins[i]

        # use a hashable container for quicker lookup
        accept_bins = set(accept_bins)

        # collect those values not in the excluded rows/columns
        keep_row = []
        keep_col = []
        keep_data = []
        for i in xrange(_map.nnz):
            if _map.row[i] in accept_bins and _map.col[i] in accept_bins:
                keep_row.append(_map.row[i])
                keep_col.append(_map.col[i])
                keep_data.append(_map.data[i])

        # after removal of those intervening, shift the coordinates of affected bins
        # TODO this could be moved into the loop above
        _shift = np.cumsum((~_order['mask']) * _bins)
        _csbins = np.cumsum(_bins)
        for i in xrange(len(keep_row)):
            # rather than build a complete list of shifts across matrix, we'll
            # sacrifice some CPU and do lookups for the appropriate bin
            ir = np.searchsorted(_csbins, keep_row[i], side='right')
            ic = np.searchsorted(_csbins, keep_col[i], side='right')
            keep_row[i] -= _shift[ir]
            keep_col[i] -= _shift[ic]

        return sp.coo_matrix((keep_data, (keep_row, keep_col)), shape=_map.shape - _shift[-1])

    def plot_seqnames(self, fname, simple=True, permute=False, **kwargs):
        """
        Plot the contact map, annotating the map with sequence names. WARNING: This can often be too dense
        to be legible when there are many (1000s) of sequences.

        :param fname: output file name
        :param simple: True plot seq map, False plot the extent map
        :param permute: permute the map with the present order
        :param kwargs: additional options passed to plot()
        """
        if permute:
            seq_id_iter = self.order.accepted_positions()
        else:
            seq_id_iter = xrange(self.order.count_accepted())

        tick_labs = []
        for i in seq_id_iter:
            if self.order.order[i]['ori'] < 0:
                tick_labs.append('- {}'.format(self.seq_info[i].name))
            else:
                tick_labs.append('+ {}'.format(self.seq_info[i].name))

        if simple:
            step = 2 if self.is_tipbased() else 1
            tick_locs = xrange(2, step*self.order.count_accepted()+step, step)
        else:
            if permute:
                _cbins = np.cumsum(self.grouping.bins[self.order.accepted_positions()])
            else:
                _cbins = np.cumsum(self.grouping.bins[self.order.accepted()])
            tick_locs = _cbins - 0.5

        self.plot(fname, permute=permute, simple=simple, tick_locs=tick_locs, tick_labs=tick_labs, **kwargs)

    def plot(self, fname, simple=False, tick_locs=None, tick_labs=None, norm=True, permute=False, pattern_only=False,
             dpi=180, width=25, height=22, zero_diag=None, alpha=0.01, robust=False, max_image_size=None,
             flatten=False):
        """
        Plot the contact map. This can either be as a sparse pattern (requiring much less memory but without visual
        cues about intensity), simple sequence or full binned map and normalized or permuted.

        :param fname: output file name
        :param tick_locs: major tick locations (minors take the midpoints)
        :param tick_labs: minor tick labels
        :param simple: if true, sequence only map plotted
        :param norm: normalize intensities by geometric mean of lengths
        :param permute: reorder map to current order
        :param pattern_only: plot only a sparse pattern (much lower memory requirements)
        :param dpi: adjust DPI of output
        :param width: plot width in inches
        :param height: plot height in inches
        :param zero_diag: set bright self-interactions to zero
        :param alpha: log intensities are log (x + alpha)
        :param robust: use seaborn robust dynamic range feature
        :param max_image_size: maximum allowable image size before rescale occurs
        :param flatten: for tip-based, flatten matrix rather than marginalise
        """

        plt.style.use('ggplot')

        fig = plt.figure()
        fig.set_figwidth(width)
        fig.set_figheight(height)
        ax = fig.add_subplot(111)

        if simple or self.bin_size is None:
            # prepare the map if not already done. This overwrites
            # any current ordering mask beyond the primary acceptance mask
            if self.processed_map is None:
                self.prepare_seq_map(norm=norm, bisto=True)
            _map = self.get_subspace(permute=permute, marginalise=False if flatten else True, flatten=flatten)
            # unless requested, zero diagonal for simple plots as its intensity tends to obscure detail
            if zero_diag is None:
                _map.setdiag(0)
            # amplify values for plotting
            _map *= 100
        else:
            _map = self.get_extent_map(norm=norm, permute=permute)

        if pattern_only:
            # sparse matrix plot, does not support pixel intensity
            if zero_diag:
                _map.setdiag(0)
            ax.spy(_map.tocsr(), markersize=5 if simple else 1)

        else:
            # a dense array plot

            # if too large, reduced it while sparse.
            if max_image_size is not None:
                full_size = _map.shape
                if np.max(full_size) > max_image_size:
                    reduce_factor = int(np.ceil(np.max(full_size) / float(max_image_size)))
                    logger.info('Full {} image reduction factor: {}'.format(full_size, reduce_factor))
                    # downsample the map
                    _map = sparse_utils.downsample(_map, reduce_factor)
                    # ticks adjusted to match
                    tick_locs = np.floor(tick_locs.astype(np.float) / reduce_factor)
                    logger.info('Map has been reduced from {} to {}'.format(full_size, _map.shape))

            _map = _map.toarray()

            if zero_diag:
                logger.debug('Removing diagonal')
                np.fill_diagonal(_map, 0)

            _map = np.log(_map + alpha)

            logger.debug('Making raster image')
            seaborn.heatmap(_map, robust=robust, square=True, linewidths=0, ax=ax, cbar=False)

        if tick_locs is not None:

            plt.tick_params(axis='both', which='both',
                            right=False, left=False, bottom=False, top=False,
                            labelright=False, labelleft=False, labelbottom=False, labeltop=False)

            if tick_labs is not None:
                min_labels = ticker.FixedFormatter(tick_labs)
                ax.tick_params(axis='y', which='minor', left=True, labelleft=True, labelsize=10)

                min_ticks = ticker.FixedLocator(tick_locs[:-1] + 0.5 * np.diff(tick_locs))

                ax.yaxis.set_minor_formatter(min_labels)
                ax.yaxis.set_minor_locator(min_ticks)

            # seaborn will not display the grid, so we make our own.
            ax.hlines(tick_locs, *ax.get_xlim(), color='grey', linewidth=0.5, linestyle='-.')
            ax.vlines(tick_locs, *ax.get_ylim(), color='grey', linewidth=0.5, linestyle='-.')

        logger.debug('Saving plot')
        fig.tight_layout()
        plt.savefig(fname, dpi=dpi)
        plt.close(fig)

