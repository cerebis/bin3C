import itertools
import numpy as np
from numba import jit, vectorize, int32, float64
from numpy import log, pi

# Likelihood Pareto2 parameters
MIN_FIELD = 2e-8
P2ALPHA = 0.122123774414444
P2LAMBDA = 13.675170758388262
P2MU = 13.973247315647466


@vectorize([float64(float64)])
def piecewise_3c(s):
    pr = MIN_FIELD
    if s < 500e3:
        # Pareto2
        pr = P2ALPHA / P2LAMBDA * (1 + (s - P2MU)/P2LAMBDA)**(- P2ALPHA - 1)
    return pr


@jit(float64(int32[:, :], float64[:, :]), nopython=True)
def poisson_lpmf2(ob, ex):
    """
    Entirely skips terms where no observational counts were recorded.

    :param ob: observed counts
    :param ex: expected counts
    :return: log likelihood
    """
    s = 0.0
    for i in xrange(ob.shape[0]):
        for j in xrange(ob.shape[1]):
            aij = ob[i, j]
            bij = ex[i, j]
            if aij == 0:
                continue
            s += aij * log(aij/bij) + bij - aij + 0.5 * log(2.0 * pi * aij)
    return -s


@jit(float64(int32[:, :], float64[:, :]), nopython=True)
def poisson_lpmf3(ob, ex):
    """
    All terms calculated.

    :param ob: observed counts
    :param ex: expected counts
    :return: log likelihood
    """
    s = 0.0
    for i in xrange(ob.shape[0]):
        for j in xrange(ob.shape[1]):
            aij = ob[i, j]
            bij = ex[i, j]
            if aij == 0:
                s += bij
            else:
                s += aij * log(aij/bij) + bij - aij + 0.5 * log(2.0 * pi * aij)
    return -s


def calc_likelihood(cm):
    """
    For a given order and ContactMap instance, calculate the log likelihood.

    :param cm: instance of ContactMap with matching identifiers
    :return: log likelihood
    """

    borders = cm.grouping.borders
    centers = cm.grouping.centers
    extent_map = cm.extent_map.tocsr().astype(np.int32)
    total_obs = cm.map_weight()

    lengths = cm.order.order['length']
    ori = cm.order.order['ori']

    log_l = 0.0
    for i, j in itertools.combinations(xrange(cm.total_seq), 2):

        # inter-contig separation defined by cumulative
        # intervening contig length.
        gap_length = cm.order.intervening(i, j)

        # contig lengths
        li = lengths[i]
        lj = lengths[j]

        # bin centers for each contig, relative to middle of each contig
        c_ik = centers[i]
        c_jl = centers[j]

        # orientation of sequences
        s_i = ori[i]
        s_j = ori[j]

        # all separations between bins, including the potential intervening distance L
        d_ij = gap_length + 0.5*(li + lj) + s_i * c_jl - s_j * c_ik.T

        # conversion to expected counts
        q_ij = total_obs * piecewise_3c(d_ij)

        # matrix element range which represents cross-terms between contig i and j
        i1, i2 = borders[i]
        j1, j2 = borders[j]

        # observed counts
        # for now this is converted to dense array as we need zeros
        n_ij = extent_map[i1:i2, j1:j2].todense()

        # log likelihood
        log_l += poisson_lpmf3(n_ij, q_ij)

    return log_l
