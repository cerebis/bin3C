#!/usr/bin/env python
import argparse
import re
import sys
from collections import defaultdict

import toml

parser = argparse.ArgumentParser('Collate Aragorn batch output for tRNA counts and write TOML to STDOUT')
parser.add_argument('-T', '--total-only', action='store_true', default=False,
                    help='Only report the total number of tRNA types',)
parser.add_argument('aragorn', nargs='?', type=argparse.FileType('r'),
                       default=sys.stdin,  help='Aragorn batch-format file')
args = parser.parse_args()

skip_pattern = re.compile(r'[0-9]+ genes found')
trna_pattern = re.compile(r'^[0-9]+\s+(tRNA-[a-zA-Z]+)\s+.*$')

trna_counts = defaultdict(int)
# with open(args.aragorn_out, 'rt') as input_h:
for line in args.aragorn:
    line = line.strip()
    if not line:
        break
    if line.startswith('>') or skip_pattern.match(line) is not None:
        continue
    m = trna_pattern.match(line)
    if m is not None:
        trna_counts[m.group(1)] += 1

if args.total_only:
    print(len(trna_counts))
else:
    trna_counts = dict(sorted(trna_counts.items(), key=lambda x: x[0]))
    print(toml.dumps({'total': len(trna_counts), 'copies': trna_counts}))
