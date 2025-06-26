#!/usr/bin/env python
import argparse
import sys

import pandas as pd

parser = argparse.ArgumentParser(description='Select contaminated bins from consolidated QC table')
parser.add_argument('--method', choices=['CheckMv1', 'CheckMv2', 'CoCoPye'], default='CoCoPye',
                    help='Method to use for contamination assessment')
parser.add_argument('--min-extent', default=500_000, type=int,
                    help='Minimum bin extent to be considered (default: 500,000 bp)')
parser.add_argument('--min-completeness', default=50, type=int,
                    help='Minimum bin completeness to be considered (default: 50)')
parser.add_argument('--min-contamination', default=10, type=int,
                    help='Minimum bin contamination to be considered (default: 10)')
parser.add_argument('QC_TABLE', help='Path to consolidated QC table')
args = parser.parse_args()

df = pd.read_csv(args.QC_TABLE, header=[0, 1], index_col=0)
df.sort_index(axis=1, inplace=True)
accepted = df.query('@df.bin3C.extent > @args.min_extent').index
df_qc = df.loc[accepted, (args.method, )] \
    .query('Completeness > @args.min_completeness and Contamination > @args.min_contamination') \
    .index \
    .to_frame() \
    .to_csv(sys.stdout, index=False, header=False)
