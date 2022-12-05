#!/usr/bin/env python

import sys
from itertools import islice

print ('\t'.join(["CHROM", "POS", "AN_ACB", "AC_ACB", "AF_ACB", "MAF_ACB", "AN_ASW", "AC_ASW", "AF_ASW", "MAF_ASW", "AN_BEB", "AC_BEB", "AF_BEB", "MAF_BEB", "AN_CDX", "AC_CDX", "AF_CDX", "MAF_CDX", "AN_CEU", "AC_CEU", "AF_CEU", "MAF_CEU", "AN_CHB.1KGP", "AC_CHB.1KGP", "AF_CHB.1KGP", "MAF_CHB.1KGP", "AN_CHN.NyuWa", "AC_CHN.NyuWa", "AF_CHN.NyuWa", "MAF_CHN.NyuWa", "AN_CHS.1KGP", "AC_CHS.1KGP", "AF_CHS.1KGP", "MAF_CHS.1KGP", "AN_CHS.NyuWa", "AC_CHS.NyuWa", "AF_CHS.NyuWa", "MAF_CHS.NyuWa", "AN_CLM", "AC_CLM", "AF_CLM", "MAF_CLM", "AN_ESN", "AC_ESN", "AF_ESN", "MAF_ESN", "AN_FIN", "AC_FIN", "AF_FIN", "MAF_FIN", "AN_GBR", "AC_GBR", "AF_GBR", "MAF_GBR", "AN_GIH", "AC_GIH", "AF_GIH", "MAF_GIH", "AN_GWD", "AC_GWD", "AF_GWD", "MAF_GWD", "AN_IBS", "AC_IBS", "AF_IBS", "MAF_IBS", "AN_ITU", "AC_ITU", "AF_ITU", "MAF_ITU", "AN_JPT", "AC_JPT", "AF_JPT", "MAF_JPT", "AN_KHV", "AC_KHV", "AF_KHV", "MAF_KHV", "AN_LWK", "AC_LWK", "AF_LWK", "MAF_LWK", "AN_MSL", "AC_MSL", "AF_MSL", "MAF_MSL", "AN_MXL", "AC_MXL", "AF_MXL", "MAF_MXL", "AN_PEL", "AC_PEL", "AF_PEL", "MAF_PEL", "AN_PJL", "AC_PJL", "AF_PJL", "MAF_PJL", "AN_PUR", "AC_PUR", "AF_PUR", "MAF_PUR", "AN_STU", "AC_STU", "AF_STU", "MAF_STU", "AN_TSI", "AC_TSI", "AF_TSI", "MAF_TSI", "AN_YRI", "AC_YRI", "AF_YRI", "MAF_YRI"]))
for fname in sys.argv[1:]:
    with open(fname, 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            print ('\t'.join(line))

