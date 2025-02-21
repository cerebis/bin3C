#!/bin/env python
import pandas as pd
import numpy as np
import Bio.SeqIO as SeqIO
from pyvenn import venn
import os
import matplotlib.pyplot as plt
import errno

from pandas.errors import EmptyDataError

try:
    from Bio.SeqUtils import gc_fraction
    def GC(sequence):
        return 100 * gc_fraction(sequence, ambiguous="ignore")
except ImportError:
    # Older versions have this:
    from Bio.SeqUtils import GC


def add_classifier(df, col_name, hit_index):
    df[col_name] = False
    df.loc[hit_index, col_name] = True


def test_existence(fn):
    if not os.path.exists(fn):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), fn)


def plot_venn(sets, names):
    assert len(sets) == len(names), 'sets and names must be the same length'

    if len(sets) == 2:
        fig, ax = venn.venn2(
            venn.get_labels(sets),
            names=names)
    elif len(sets) == 3:
        fig, ax = venn.venn3(
            venn.get_labels(sets),
            names=names)
    elif len(sets) == 4:
        fig, ax = venn.venn4(
            venn.get_labels(sets),
            names=names)
    elif len(sets) == 5:
        fig, ax = venn.venn5(
            venn.get_labels(sets),
            names=names)
    else:
        raise RuntimeError("number of sets must be between 2 and 5")
    return fig, ax


def combine_results(out_dir, contig_fasta, circ_file, cat_file,
                    virsort_file, vibrant_file, dvf_file, marvel_file, phamer_file,
                    plasflow_file, plasclass_file, plasforest_file, plasmidhunter_file,
                    no_partial=True, silent=False):

    # check for input and output existence
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        print(f'Created output directory: {out_dir}')

    for fn in [contig_fasta, virsort_file, vibrant_file, dvf_file, marvel_file, phamer_file,
               plasflow_file, plasclass_file, plasforest_file]:
        assert os.path.exists(fn), f'Input file {fn} does not exist'

    base_name = os.path.splitext(os.path.basename(contig_fasta))[0]

    if not silent:
        print('Parsing contigs for sequence details...')
    all_seq = []
    for _s in SeqIO.parse(contig_fasta, 'fasta'):
        all_seq.append({'name': _s.id, 'seq_len': len(_s), 'gc': GC(_s.seq)})
    all_seq = pd.DataFrame(all_seq).set_index('name')
    all_seq.insert(loc=0, column='run', value=base_name)
    if not silent:
        print('There are {} contigs in the assembly'.format(len(all_seq)))

    if cat_file is not None:
        cat = pd.read_csv(cat_file, sep='\t')
        cat.rename(columns={'# contig': 'name'}, inplace=True)
        cat.set_index('name', inplace=True)
        cat[['domain','phylum','class','order','family','genus','species']] = cat.lineage.str.split(';', expand=True).iloc[:,1:]
        #        cat['superkingdom.1'] = cat.superkingdom.str.split(':', expand=True)[0]
        #        cat['phylum.1'] = cat.phylum.str.split(':', expand=True)[0]
        all_seq = all_seq.join(cat, how='left')
        print('{:,} contig annotations'.format(len(cat)))

    if circ_file is not None:
        circ = set(pd.read_csv(circ_file, index_col=0, names=['seq']).index)
        ix_circ = all_seq.index.isin(circ)
        all_seq['gfa_topology'] =  'linear'
        all_seq.loc[ix_circ, 'gfa_topology'] = 'circular'

    # VIRUS results

    virsort = pd.read_csv(virsort_file, sep='\t')
    split_col = virsort.seqname.str.extract('^([^\|]+).*?([^\|]+)')
    virsort['name'] = split_col[0]
    virsort['virsort_partial'] = ~ split_col[1].str.startswith('full')
    # confidence score, or a hallmark and viral>cellular content and only full length
    virsort.query('(max_score>0.7 or (hallmark > 0 and viral > cellular)) and length>=5000', inplace=True)
    if no_partial:
        virsort.query('~ virsort_partial', inplace=True)
    virsort.set_index('name', inplace=True)
    add_classifier(all_seq, 'virsorter', virsort.index)
    # record virsort group type
    all_seq['virsort_group'] = None
    all_seq.loc[virsort.index, 'virsort_group'] = virsort['max_score_group']
    # record virsort partial/full determination
    all_seq['virsort_partial'] = None
    all_seq.loc[virsort.index, 'virsort_partial'] = virsort['virsort_partial']
    if not silent:
        print('{:,} virsorter predictions'.format(len(virsort)))

    try:
        vib = pd.read_csv(vibrant_file, header=None)
        vib = vib[0].str.split('_fragment_', expand=True)
        vib[0] = vib[0].str.split(' ', expand=True)[0]
        vib = vib.rename(columns={0: 'name'})
        if len(vib.columns) == 1:
            # handle case when no seq was marked as a fragment
            vib['vibrant_partial'] = None
        else:
            # set partial flag
            vib = vib.rename(columns={1: 'vibrant_partial'})
            vib = vib.assign(vibrant_partial = lambda x: ~pd.isna(x.vibrant_partial))
        vib.set_index('name', inplace=True)

        if not silent:
            print('vibrant: {} initial hits'.format(len(vib)))
        if no_partial:
            n_frag = vib.vibrant_partial.sum()
            vib = vib.query('~ vibrant_partial')
            if not silent:
                print('vibrant: dropped {} fragments'.format(n_frag))
    except EmptyDataError:
        # handle case when vibrant has to predictions (0-byte file)
        vib = pd.DataFrame({'name': [], 'vibrant_partial': []}, dtype=str).set_index('name')

    add_classifier(all_seq, 'vibrant', vib.index)
    all_seq['vibrant_partial'] = None
    all_seq.loc[vib.index, 'vibrant_partial'] = vib['vibrant_partial']
    if not silent:
        print('{:,} vibrant predictions'.format(len(vib)))

    dvf = pd.read_csv(dvf_file, sep='\t')
    dvf = dvf.query('score >= 0.85 and pvalue <= 0.05')
    dvf.set_index('name', inplace=True)
    add_classifier(all_seq, 'deepvirfinder', dvf.index)
    if not silent:
        print('{:,} deepvirfinder predictions'.format(len(dvf)))

    marv = pd.read_csv(marvel_file, sep='\t', header=None) \
        .rename(columns={0: 'name', 1: 'PHAGE_PROB', 2: 'IS_PHAGE'}) \
        .set_index('name')
    marv = marv.query('PHAGE_PROB >= 0.8')
    add_classifier(all_seq, 'marvel', marv.index)
    if not silent:
        print('{:,} marvel predictions'.format(len(marv)))

    phamer = pd.read_csv(phamer_file, sep='\t').query('Pred == "virus" and PhaMerConfidence == "high-confidence"')
    phamer.set_index('Accession', inplace=True)
    add_classifier(all_seq, 'phamer', phamer.index)
    if not silent:
        print('{:,} phamer predictions'.format(len(phamer)))

    # create sets of contig ids for each 2-way consensus
    ix1 = set(virsort.join(dvf, how='inner').index)
    ix2 = set(virsort.join(vib, how='inner').index)
    ix3 = set(virsort.join(marv, how='inner').index)
    ix4 = set(virsort.join(phamer, how='inner').index)
    ix5 = set(dvf.join(vib, how='inner').index)
    ix6 = set(dvf.join(marv, how='inner').index)
    ix7 = set(dvf.join(phamer, how='inner').index)
    ix8 = set(vib.join(marv, how='inner').index)
    ix9 = set(vib.join(phamer, how='inner').index)
    ix10 = set(marv.join(phamer, how='inner').index)

    # matrix representation
    if not silent:
        df_2way = np.zeros((5,5), dtype=np.int32)
        df_2way[np.triu_indices(5,k=1)] = [len(ix1),len(ix2),len(ix3),len(ix4),len(ix5),len(ix6),len(ix7),len(ix8),len(ix9),len(ix10)]
        df_2way[np.diag_indices(5)] = [len(virsort), len(dvf), len(vib), len(marv), len(phamer)]
        df_2way = pd.DataFrame(df_2way, columns=['vs2', 'dvf', 'vib', 'marv', 'phamer'], index=['vs2', 'dvf', 'vib', 'marv', 'phamer'])
        print('Two-way consensus matrix')
        print(df_2way)


    ax = plt.subplot(111)
    fig, ax = plot_venn([set(virsort.index), set(dvf.index), set(marv.index), set(vib.index), set(phamer.index)],
                        ['virsorter', 'deepvirfinder', 'marvel', 'vibrant', 'phamer'])
    fig = ax.figure
    fig.set_figwidth(10)
    fig.set_figheight(10)
    fig.savefig(os.path.join(out_dir, f'{base_name}_virus-venn.svg'), bbox_inches='tight')
    plt.close()

    all_seq['virus_hit_rank'] = all_seq.loc[:, ['virsorter', 'deepvirfinder', 'vibrant', 'marvel', 'phamer']].sum(axis=1)
    vir_2way = all_seq.query('virus_hit_rank>=2')
    if not silent:
        print('{:,} 2-way consensus virus identifications\n'.format(len(vir_2way)))

    vir_2way.to_csv(os.path.join(out_dir, f'{base_name}_virus_2way.csv'))
    vir_2way.reset_index().name.to_csv(os.path.join(out_dir, f'{base_name}_virus_2way.ids'), index=False, header=False)

    # PLASMID results

    plasflow = pd.read_csv(plasflow_file, sep='\t', na_values='None').rename(columns={'contig': 'name'})
    plasflow.query('score>=0.9 and classification.str.startswith("plasmid")', engine='python', inplace=True)
    plasflow.set_index('name', inplace=True)
    add_classifier(all_seq, 'plasflow', plasflow.index)
    if not silent:
        print('{:,} plasflow predictions'.format(len(plasflow)))

    plasclass = pd.read_csv(plasclass_file, sep='\t', names=['contig_id', 'pc_prob'])
    plasclass.query('pc_prob>=0.9', inplace=True)
    plasclass.set_index('contig_id', inplace=True)
    add_classifier(all_seq, 'plasclass', plasclass.index)
    if not silent:
        print('{:,} plasclass predictions'.format(len(plasclass)))

    plasforest = pd.read_csv(plasforest_file, sep=',').rename(columns={'ID': 'name'})
    plasforest.query('Prediction == "Plasmid"', inplace=True)
    plasforest.set_index('name', inplace=True)
    add_classifier(all_seq, 'plasforest', plasforest.index)
    if not silent:
        print('{:,} plasforest predictions'.format(len(plasforest)))

    plashunt = pd.read_csv(plasmidhunter_file, sep='\t').rename(columns={'Unnamed: 0': 'seq',
                                                                         'Prediction (0: chromosome, 1: plasmid)': 'Pred'})
    plashunt.set_index('seq', inplace=True)
    plashunt.query('Pred == 1', inplace=True)
    add_classifier(all_seq, 'plasmidhunter', plashunt.index)
    if not silent:
        print('{:,} plasmidhunter predictions'.format(len(plashunt)))

    # create sets of contig ids for each 2-way consensus
    ix1 = set(plasflow.join(plasclass, how='inner').index)
    ix2 = set(plasflow.join(plasforest, how='inner').index)
    ix3 = set(plasflow.join(plashunt, how='inner').index)
    ix4 = set(plasclass.join(plasforest, how='inner').index)
    ix5 = set(plasclass.join(plashunt, how='inner').index)
    ix6 = set(plasforest.join(plashunt, how='inner').index)

    # matrix representation
    if not silent:
        df_2way = np.zeros((4,4), dtype=np.int32)
        df_2way[np.triu_indices(4,k=1)] = [len(ix1),len(ix2),len(ix3),len(ix4),len(ix5),len(ix6)]
        df_2way[np.diag_indices(4)] = [len(plasflow), len(plasclass), len(plasforest), len(plashunt)]
        df_2way = pd.DataFrame(df_2way, columns=['pflw', 'pcls', 'pfor', 'phnt'], index=['pflw', 'pcls', 'pfor', 'phnt'])
        print('Two-way consensus matrix')
        print(df_2way)

    ax = plt.subplot(111)
    fig, ax = plot_venn([set(plasflow.index), set(plasclass.index), set(plasforest.index), set(plashunt.index)],
                        ['plasflow', 'plasclass', 'plasforest', 'plasmidhunter'])
    fig = ax.figure
    fig.set_figwidth(10)
    fig.set_figheight(10)
    fig.savefig(os.path.join(out_dir, f'{base_name}_plasmid-venn.svg'), bbox_inches='tight')
    plt.close()

    all_seq['plasmid_hit_rank'] = all_seq.loc[:, ['plasflow', 'plasclass', 'plasforest', 'plasmidhunter']].sum(axis=1)
    pls_2way = all_seq.query('plasmid_hit_rank>=2')
    if not silent:
        print('{:,} 2-way consensus plasmid identifications'.format(len(pls_2way)))

    pls_2way.to_csv(os.path.join(out_dir, f'{base_name}_plasmid_2way.csv'))
    pls_2way.reset_index().name.to_csv(os.path.join(out_dir, f'{base_name}_plasmid_2way.ids'), index=False, header=False)

    # OVERALL results

    all_seq['mge_status'] = 'chromosome'
    all_seq.loc[(all_seq.virus_hit_rank>=2) & (all_seq.plasmid_hit_rank<2), 'mge_status'] = 'virus'
    all_seq.loc[(all_seq.virus_hit_rank<2) & (all_seq.plasmid_hit_rank>=2), 'mge_status'] = 'plasmid'
    all_seq.loc[(all_seq.virus_hit_rank>=2) & (all_seq.plasmid_hit_rank>=2), 'mge_status'] = 'mge'

    # venn showing cross-over of hits between viruses and plasmids
    ix_v2w = set(all_seq.query('virus_hit_rank>=2').index)
    ix_p2w = set(all_seq.query('plasmid_hit_rank>=2').index)
    ax = plt.subplot(111)
    fig, ax = plot_venn([ix_v2w, ix_p2w], ['Virus 2-way', 'Plasmid 2-way'])
    fig = ax.figure
    fig.set_figwidth(10)
    fig.set_figheight(10)
    fig.savefig(os.path.join(out_dir, f'{base_name}_virus_to_plasmid-venn.svg'), bbox_inches='tight')
    plt.close()
    if not silent:
        print('\n{:,} sequences were 2-way consensus hits for both virus and plasmid'.format(
            len(all_seq.query('virus_hit_rank>=2 and plasmid_hit_rank>=2'))))

    mge_2way = all_seq[all_seq.index.isin(ix_v2w | ix_p2w)]
    mge_2way.to_csv(os.path.join(out_dir, f'{base_name}_mge_2way.csv'))
    mge_2way.reset_index().name.to_csv(os.path.join(out_dir, f'{base_name}_mge_2way.ids'), index=False, header=False)

    all_seq.to_csv(os.path.join(out_dir, f'{base_name}_report.csv'))

    return all_seq

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Collate MGE predictions from multiple tools')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('contig_fasta', help='Contig fasta file')
    parser.add_argument('circ_file', help='Circular contig file')
    parser.add_argument('cat_file', help='Contig annotation file')
    parser.add_argument('virsort_file', help='VirSorter output file')
    parser.add_argument('vibrant_file', help='VIBRANT output file')
    parser.add_argument('dvf_file', help='DeepVirFinder output file')
    parser.add_argument('marvel_file', help='MARVEL output file')
    parser.add_argument('phamer_file', help='PhaMer output file')
    parser.add_argument('plasflow_file', help='PlasFlow output file')
    parser.add_argument('plasclass_file', help='PlasClass output file')
    parser.add_argument('plasforest_file', help='PlasForest output file')
    parser.add_argument('plasmidhunter_file', help='PlasmidHunter output file')
    parser.add_argument('--no_partial', action='store_true', help='Ignore partial predictions')
    parser.add_argument('--silent', action='store_true', help='Suppress output')

    args = parser.parse_args()
    combine_results(args.out_dir, args.contig_fasta, args.circ_file, args.cat_file,
                    args.virsort_file, args.vibrant_file, args.dvf_file, args.marvel_file, args.phamer_file,
                    args.plasflow_file, args.plasclass_file, args.plasforest_file, args.plasmidhunter_file,
                    no_partial=args.no_partial, silent=args.silent)
