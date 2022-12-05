#!/usr/bin/env python

import sys
from itertools import islice

print ('\t'.join(["CHROM", "POS", "AN_AFR", "AC_AFR", "AF_AFR", "MAF_AFR", "AN_AMR", "AC_AMR", "AF_AMR", "MAF_AMR", "AN_EAS", "AC_EAS", "AF_EAS", "MAF_EAS", "AN_EUR", "AC_EUR", "AF_EUR", "MAF_EUR", "AN_SAS", "AC_SAS", "AF_SAS", "MAF_SAS"]))
for fname in sys.argv[1:]:
    with open(fname, 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            print ('\t'.join(line))

