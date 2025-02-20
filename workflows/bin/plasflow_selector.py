#!/bin/env python
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--threshold', type=float, help='Probability threshold', default=0.9)
parser.add_argument('score_file', help='PlasFlow score output file')
args = parser.parse_args()

assert 0 <= args.threshold <= 1, 'Probability threshold must be in the range [0, 1]'

df = pd.read_csv(args.score_file, sep='\t')
df.drop(columns=df.columns[0], inplace=True)
df.set_index('contig_id', inplace=True)

print('contig\tclassification\tscore')

for index, row in df.iterrows():
    label_name = pd.to_numeric(row[4:]).idxmax()
    tax_name = label_name.split(".", 1)[1]

    if row[label_name] < args.threshold:

        plasmids_sum = row[[col for col in df.columns if re.match(r'^plasmid.*', col)]].sum()
        chromosomes_sum = row[[col for col in df.columns if re.match(r'^chromosom.*', col)]].sum()

        my_regex = r".*" + re.escape(tax_name) + r""
        taxnames_sum = row[[col for col in df.columns if re.match(my_regex, col)]].sum()

        if plasmids_sum > args.threshold:
            assignment = 'plasmid.unclassified'
            score = plasmids_sum
        elif chromosomes_sum > args.threshold:
            assignment =  'chromosome.unclassified'
            score = chromosomes_sum
        elif taxnames_sum > args.threshold:
            assignment =  'unclassified.{}'.format(tax_name)
            score = taxnames_sum
        else:
            assignment =  'unclassified.unclassified'
            score = None
    else:
        assignment = label_name
        score = row[label_name]

    print('{}\t{}\t{}'.format(row['contig_name'], assignment, score))
