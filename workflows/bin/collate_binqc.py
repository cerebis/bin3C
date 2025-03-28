#!/usr/bin/env python
import os
import pandas as pd
import re
from collections import defaultdict
from plotnine import *
from pypalettes import load_cmap


# point of brittle failure
QC_METHODS = ['CheckMv1', 'CheckMv2', 'CoCoPye']


def stat_larger(df, min_completeness=50, max_contamination=500):
    return df.query('Completeness >= @min_completeness and Contamination <= @max_contamination')

def read_gtdb(fname, drop_dupes=False):
    df = pd.read_csv(fname, sep='\t') \
        .rename(columns={'# bin': 'Name', 'lineage scores (f: 0.30)': 'lineage_scores'})
    if drop_dupes:
        df = df.drop_duplicates(subset='Name', keep='first')
    df['Name'] = df.Name.str.split('.', expand=True)[0]
    return df

def count_tRNA(aragorn_file):
    """
    Count tRNA genes in an Aragorn output file.

    Tabulate the number of copies of each tRNA gene type. Return the total number of
    types detected by Aragorn.
    :param aragorn_file:
    :return: dict reporting the total number of tRNA gene types (tRNA-[3let]) and the number of copies of each type
    """
    skip_pattern = re.compile(r'[0-9]+ genes found')
    trna_pattern = re.compile(r'^[0-9]+\s+(tRNA-[a-zA-Z]+)\s+.*$')
    trna_counts = defaultdict(int)
    with open(aragorn_file, 'rt') as input_h:
        for line in input_h:
            line = line.strip()
            if not line:
                break
            if line.startswith('>') or skip_pattern.match(line) is not None:
                continue
            m = trna_pattern.match(line)
            trna_counts[m.group(1)] += 1
    trna_counts = dict(sorted(trna_counts.items(), key=lambda x: x[0]))
    return {'total': len(trna_counts), 'copies': trna_counts}


def count_rRNA(barrnap_file):
    """
    Count rRNA genes in a Barrnap output file.

    The expectation is that we are analysing prokaryotes and
    therefore only 5S, 16S and 23S are of interest, however Archaeal and Eukaryotic rRNAs are also included.

    :param barrnap_file:
    :return: Series with counts of each type of rRNA gene
    """
    df = pd.read_csv(barrnap_file, sep='\t', comment='#', header=None)
    df = df.query('not @df[8].str.contains("partial")')[8].str.extract(r'Name=(\w+)')
    df.columns = ['gene']
    return pd.Categorical(df.gene, categories=['5S_rRNA', '16S_rRNA', '23S_rRNA',
                                               '5_8S_rRNA', '12S_rRNA', '18S_rRNA', '28S_rRNA']).value_counts()

def assign_mimag_qualities(df):

    def select_quality_label(row, method):

        # forced to sort column index each row to suppress a warning.
        #   sorting the parent table seems to have no impact despite too
        #   much time spent troubleshooting.
        row = row.sort_index()

        # QC values from the given method
        qc = row[(method,)]
        # RNA results
        rna = row[('RNA_report',)]

        # anything over 10% is a failure
        # CoCoPye also assigns -100 when it abandons a bin
        if qc.Contamination > 10 or (qc.Contamination < 0 or qc.Completeness < 0):
            return 'Failed'

        if qc.Completeness >= 90 and qc.Contamination <= 5:
            if rna.tRNA_count >= 18 and rna['5S_rRNA']>0 and rna['16S_rRNA']>0 and rna['23S_rRNA']>0:
                return 'High'
            else:
                # an addition non-standard label, when not all RNA are present
                return 'Near'
        elif qc.Completeness >= 50 and qc.Contamination <= 10:
            return 'Medium'
        elif qc.Completeness < 50 and qc.Contamination <= 10:
            return 'Low'
        # ideally there should be no None
        return '(not assigned)'

    cats = ['Failed','Low','Medium','Near','High']
    for _method in QC_METHODS:
        quals = df.apply(select_quality_label, axis=1, method=_method)
        df[(_method, 'MIMAG_quality')] = pd.Categorical(quals, categories=cats)

def combine_qc_results(out_dir,
                       cluster_report, rna_report,
                       checkm1_file, checkm2_file, cocopye_file,
                       gtdbtk_file):

    def add_fancy_columns(df, name):
        df.columns = pd.MultiIndex.from_product([[name], df.columns])

    df_report = pd.read_csv(cluster_report) \
        .rename(columns={'name': 'bin'}) \
        .set_index('bin')
    add_fancy_columns(df_report, 'bin3C')

    df_rna = pd.read_csv(rna_report) \
        .set_index('bin')
    add_fancy_columns(df_rna, 'RNA_report')

    df_chkm1 = pd.read_csv(checkm1_file, sep='\t') \
        .rename(columns={'Bin Id': 'bin'}) \
        .set_index('bin')
    add_fancy_columns(df_chkm1, 'CheckMv1')

    df_chkm2 = pd.read_csv(checkm2_file, sep='\t') \
        .rename(columns={'Name': 'bin'}) \
        .set_index('bin')
    add_fancy_columns(df_chkm2, 'CheckMv2')

    df_coco = pd.read_csv(cocopye_file) \
        .assign(completeness=lambda x: x.completeness * 100,
                contamination=lambda x: x.contamination * 100) \
        .rename(columns={'completeness': 'Completeness',
                         'contamination': 'Contamination'}) \
        .set_index('bin')
    add_fancy_columns(df_coco, 'CoCoPye')

    df_gtdb = pd.read_csv(gtdbtk_file, sep='\t') \
        .rename(columns={'user_genome': 'bin'}) \
        .drop_duplicates(subset='bin', keep='first') \
        .set_index('bin')
    df_gtdb[['domain','phylum','class','order','family','genus','species']] = df_gtdb.classification.str.split(';', expand=True) \
        .apply(lambda x: [xi[3:] if xi is not None else '-' for xi in x if xi != 'root'], axis=1, result_type='expand')
    add_fancy_columns(df_gtdb, 'GTDBtk')

    # TODO refer to figures. Do we just concat instead.
    df_all = df_report \
        .join(df_rna, how='left', validate='1:1') \
        .join(df_chkm1, how='left', validate='1:1') \
        .join(df_chkm2, how='left', validate='1:1') \
        .join(df_coco, how='left', validate='1:1') \
        .join(df_gtdb, how='left', validate='1:1') #\

    df_all.sort_index()
    assign_mimag_qualities(df_all)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        print(f'Created output directory: {out_dir}')

    df_all.to_csv(f'{out_dir}/qc_summary.csv')

    # plot
    df_all.sort_index(axis=1, inplace=True, level=0)
    df_big = df_all.query('@df_all.bin3C.extent>100_000')
    qual_table = []
    for _m in QC_METHODS:
        a = df_big[(_m, )].loc[:, ['MIMAG_quality']]
        a['method'] = _m
        qual_table.append(a)
    qual_table = pd.concat(qual_table)
    qual_table.reset_index(inplace=True)

#, position=position_dodge(preserve='single'))
    p = ( ggplot(qual_table) + geom_bar(aes(x='method', fill='MIMAG_quality'), width=0.67)
          + scale_fill_manual(['#858585FF'] + load_cmap('miami2').colors[1:])
          + labs(x='QC Method', y='Number of MAGs', fill='MIMAG Quality')
          + theme_bw() + theme(figure_size=[6,5], legend_position='bottom', legend_title=element_blank())
          )
    p.save(filename=f'{out_dir}/quality_breakdown.png', dpi=300, verbose=False)
    p.save(filename=f'{out_dir}/quality_breakdown.svg', verbose=False)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser('Collate Genome Binning QC results')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('cluster_report', help='Cluster report')
    parser.add_argument('rna_report', help='RNA report')
    parser.add_argument('checkm1_file', help='CheckMv1 report')
    parser.add_argument('checkm2_file', help='CheckMv2 report')
    parser.add_argument('cocopye_file', help='CoCoPye report')
    parser.add_argument('gtdbtk_file', help='GTDBtk report')
    args = parser.parse_args()

    combine_qc_results(args.out_dir, args.cluster_report, args.rna_report,
                       args.checkm1_file, args.checkm2_file, args.cocopye_file,
                       args.gtdbtk_file)
