#!/usr/bin/env python
import json
import os
import sys

scale_factor = float(sys.argv[1])
assert scale_factor >= 1, 'Scale factor must be greater than or equal to 1'

qc3c_filename = sys.argv[2]
assert os.path.exists(qc3c_filename), f'The file {qc3c_filename} does not exist'

with open(qc3c_filename, 'rt') as input_handle:
    qc3c_report = json.load(input_handle)
    print(f'{scale_factor * qc3c_report["obs_insert_mean"]:.0f}')
