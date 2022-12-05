#!/usr/bin/env python

import sys
from itertools import islice

print ('\t'.join(["CHROM", "POS", "AN", "AC", "AF", "MAF", "AN_1KGP", "AC_1KGP", "AF_1KGP", "MAF_1KGP", "AN_NyuWa", "AC_NyuWa", "AF_NyuWa", "MAF_NyuWa"]))
for fname in sys.argv[1:]:
    with open(fname, 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            print ('\t'.join(line))

