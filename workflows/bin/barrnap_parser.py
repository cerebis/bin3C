#!/usr/bin/env python
import argparse
import sys

import pandas as pd

parser = argparse.ArgumentParser('Collate Aragorn batch output for tRNA counts and write TOML to STDOUT')
parser.add_argument('--name', type=str, help='Name of analysed sequence')
parser.add_argument('barrnap_gff', type=argparse.FileType('r'), nargs='?', default=sys.stdin,
                    help='Barrnap GFF output file')
args = parser.parse_args()

try:
    df = pd.read_csv(args.barrnap_gff, sep='\t', comment='#', header=None)
    df = df.query('not @df[8].str.contains("partial")')[8].str.extract(r'Name=(\w+)')
    df.columns = ['gene']
    out = pd.Categorical(df.gene, categories=['5S_rRNA', '16S_rRNA', '23S_rRNA',
                                               '5_8S_rRNA', '12S_rRNA', '18S_rRNA', '28S_rRNA']).value_counts()
except pd.errors.EmptyDataError:
    out = pd.Categorical([], categories=['5S_rRNA', '16S_rRNA', '23S_rRNA',
                                         '5_8S_rRNA', '12S_rRNA', '18S_rRNA', '28S_rRNA']).value_counts()

out = out.to_frame().T
if args.name:
    out['name'] = args.name

out.to_csv(sys.stdout, sep=',', index=False, header=True)
