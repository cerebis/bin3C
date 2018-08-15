#!/usr/bin/env python
import community as com
import logging
import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)


def decompose_graph(g):
    """
    Using the Louvain algorithm for community detection, as
    implemented in the community module, determine the partitioning
    which maximises the modularity. For each individual partition
    create the sub-graph of g

    :param g: the graph to decompose
    :return: the set of sub-graphs which form the best partitioning of g
    """
    decomposed = []
    p = com.best_partition(g)
    partition_labels = np.unique(p.values())

    # for each partition, create the sub-graph
    for pi in partition_labels:
        # start with a complete copy of the graph
        gi = g.copy()
        # build the list of nodes not in this partition and remove them
        to_remove = [n for n in g.nodes_iter() if p[n] != pi]
        gi.remove_nodes_from(to_remove)
        decomposed.append(gi)

    return decomposed


def cluster(g, no_iso, method=None, ragbag=False, verbose=False):

    assert not (no_iso and ragbag), 'options no_iso and ragbag are mutually exclusive'

    ragbag_group = None
    singletons = None

    # if we need them, determine singleton nodes, but don't count self-loops (which networkx does)
    if no_iso or ragbag:
        # this is a memory costly approach.
        g_nsl = g.copy()
        g_nsl.remove_edges_from(g_nsl.selfloop_edges())
        singletons = nx.isolates(g_nsl)
        g_nsl.clear()

    # remove isolated nodes and forget about em
    if no_iso:
        logger.info('Removed {0} isolated nodes from graph'.format(len(singletons)))
        g.remove_nodes_from(singletons)
        print_info(g)

    # put them in a ragbag instead
    elif ragbag:
        if len(singletons) == 0:
            logger.info('Ragbag cluster would be empty, so excluded')
            ragbag_group = {}
        else:
            logger.info('Ragbag cluster will contain {0} nodes'.format(len(singletons)))
            g.remove_nodes_from(singletons)
            ragbag_group = dict((n, 1.0) for n in singletons)
            print_info(g)

    # determine the best partitioning of g
    logger.info('Determining best partitioning')
    partitions = com.best_partition(g)

    # build a dictionary of classification from the partitioning result
    # this is effectively a hard clustering answer to the problem
    com_ids = set(partitions.values())
    logger.info('There were {0} communities in decomposition'.format(len(com_ids)))

    logger.info('Inverting partition map')
    revpart = {}
    for ni, ci in partitions.iteritems():
        revpart.setdefault(ci, []).append((ni, 1.0))

    if verbose:
        logger.info('Communities with more than one node:')

    # dict of communities
    communities = {}
    for ci in com_ids:
        communities[ci] = dict(revpart[ci])
        if verbose and len(communities[ci]) > 1:
            logger.info('\tcommunity {0}: {1} nodes'.format(ci, len(communities[ci])))

    # bit o memory conservation
    revpart.clear()

    if method == 'maxaff':
        for u in g.nodes_iter():
            if g.degree(u) > 0:
                max_u = max([x[1]['weight'] for x in g[u].items()])
                for v in nx.all_neighbors(g, u):
                    if partitions[u] != partitions[v]:
                        max_v = max([x[1]['weight'] for x in g[v].items()])
                        w_v = g[u][v]['weight']
                        if w_v == max_u:
                            communities[partitions[v]][u] = 0.5
                        if w_v == max_v:
                            communities[partitions[u]][v] = 0.5

    elif method == 'simple':
        # iterate over the nodes in the graph, if an edge connects nodes in disjoint
        # partitions, then make both nodes members of both partitions.
        for n1 in g.nodes_iter():
            for n2 in nx.all_neighbors(g, n1):
                if partitions[n1] != partitions[n2]:
                    # we add them at half weight just for discrimination
                    communities[partitions[n1]][n2] = 0.5
                    communities[partitions[n2]][n1] = 0.5

    if ragbag and len(ragbag_group) > 0:
        rb_id = max(communities)+1
        if rb_id in communities.keys():
            raise RuntimeError('ragbag id clashed. This function is really dumb and some refactoring is in order')
        communities[rb_id] = ragbag_group

    return communities


def print_info(g):
    """
    logger.info(order and size of graph g)
    :param g: graph to logger.info(info)
    """
    logger.info('Graph composed of {0} nodes and {1} edges'.format(g.order(), g.size()))


def write_mcl(communities, path):
    """
    Write communities dictionary to output file.
    :param communities: the dictionary of communities
    :param path: the output file path
    """
    with open(path, 'w') as hout:
        keys = sorted(communities.keys())
        for k in keys:
            line = ' '.join([str(sid) for sid in sorted(communities[k].keys())])
            hout.write(line.strip())
            hout.write('\n')


def write_output(communities, filename, ofmt='mcl'):
    if ofmt == 'mcl':
        write_mcl(communities, filename)
    elif ofmt == 'graphml':
        # convert communities to a graph
        cg = nx.DiGraph()
        for k, v in communities.iteritems():
            cg.add_node(k)
            for vi in v:
                cg.add_edge(k, vi)
        nx.write_graphml(cg, filename)
    else:
        raise RuntimeError('Unsupported format type: {0}'.format(ofmt))


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Decompose a graph into its communities')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Verbose output')
    parser.add_argument('--no-isolates', action='store_true', default=False,
                        help='Remove isolated nodes')
    parser.add_argument('--otype', choices=['hard', 'soft', 'maxaff'], default='hard',
                        help='Output type')
    parser.add_argument('--ifmt', choices=['edgelist', 'graphml'], default='graphml',
                        help='Specify input format [graphml]')
    parser.add_argument('--ofmt', choices=['mcl', 'graphml'], default='mcl',
                        help='Specify output format [mcl]')
    parser.add_argument('--ragbag', action='store_true', default=False,
                        help='Place isolates in a single ragbag cluster')
    parser.add_argument('input', help='Input graph (graphml format)')
    parser.add_argument('output', help='Output file')
    args = parser.parse_args()

    if args.otype == 'induced':
        raise RuntimeError('induced option no longer supported')

    if args.ifmt == 'graphml':
        g = nx.read_graphml(args.input)
    else:
        g = nx.read_edgelist(args.input, data=(('weight', float), ))

    print 'Initial statistics'
    print_info(g)

    if args.otype == 'soft':
        method = 'simple'
    elif args.otype == 'maxaff':
        method = 'maxaff'
    else:
        method = None

    communities = cluster(g, args.no_isolates, method=method, ragbag=args.ragbag, verbose=args.verbose)

    write_output(communities, args.output, args.ofmt)
