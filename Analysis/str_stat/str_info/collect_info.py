#!/usr/bin/env python

import sys
from itertools import islice


print (','.join(["CHROM", "POS", "END", "RU", "PERIOD", "REF", "EXPTHRESH", "HRUN"]))
for fname in sys.argv[1:]:
    with open(fname, 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            print (','.join(line))
