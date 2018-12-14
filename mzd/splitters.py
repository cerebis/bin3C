import numpy as np
import Bio.SeqIO as SeqIO
import logging
import os

# package logger
logger = logging.getLogger(__name__)


def simple_splitter(win_size, seq_file, out_file=None, threshold=0.333, fmt='fasta', clobber=False):
    """
    A simplistic routing which will split a set of sequences into pieces.

    Sequence IDs of the fragments will contain their relevant coordinates appended to the
    originating sequence's identifier.

    The original sequence must be long enough that the resulting fragments
    are close to the requested target length.

    :param win_size: the target size for chunks
    :param seq_file: the file name for the sequences to split
    :param out_file: if set, resulting fragmented sequences wrriten to this file name, otherwise
    the file name will be the same as input with the added suffix ".split"
    :param threshold: a wiggle factor for "nearly" long enough sequences.
    :param fmt: the format of input and output sequences.
    :param clobber: overwrite existing output files
    :return: the file name of the split sequences
    """
    if out_file is None:
        out_file = '{}.split'.format(seq_file)

    if not clobber and os.path.exists(out_file):
        raise IOError('output path already exists!')

    with open(out_file, 'w') as out_h:

        n_seqs = 0
        sum_seqs = 0
        max_seq = -1
        n_chunks = 0
        sum_x = 0
        max_x = -1
        n_x = 0

        # split longer sequences into uniformly sized pieces.
        for seq in SeqIO.parse(seq_file, fmt):

            # splitting is governed by a minimum size and threshold
            # this controls how much wiggle room is permitted for sequences
            # whose lengths are "almost" long enough for the next integer
            # number of pieces. This way, we distribute the differences
            # across pieces rather than tacking on the excess
            # onto one or two.

            l = len(seq)
            sum_seqs += l
            if l > max_seq:
                max_seq = l

            n = int(l / win_size + threshold)
            if n == 0:
                n = 1

            # determine uniformly spaced positions
            x = np.linspace(0, l, n+1, dtype=int)

            out_seqs = []
            for i in range(1, len(x)):
                s = seq[x[i-1]: x[i]]
                # output pieces include their coordinates in their identifiers
                s.id = '{}.{}_{}'.format(seq.id, x[i-1], x[i])
                out_seqs.append(s)

            dx = np.diff(x)
            sum_x += np.sum(dx)
            n_x += len(dx)
            max_dx = np.max(dx)
            if max_dx > max_x:
                max_x = max_dx

            n_chunks += len(out_seqs)
            n_seqs += 1

            # write out per input sequence
            SeqIO.write(out_seqs, out_h, fmt)

        logger.debug('There were {} input sequences with mean size {:.0f} bp and max {} bp'.format(n_seqs, sum_seqs / float(n_seqs), max_seq))
        logger.debug('Splitting produced {} fragments of mean size {:.0f} bp and max {} bp'.format(n_chunks, sum_x / float(n_x), max_x))
        logger.info('Resulting fragmented sequences written to {}'.format(out_file))

    return out_file
