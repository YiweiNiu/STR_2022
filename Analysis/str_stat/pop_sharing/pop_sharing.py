#!/usr/bin/env python

'''
purpose:
get sharing stat of pSTRs

usage:
cd /home2/niuyw/project/STR/str_stat/pop_sharing
python3 pop_sharing.py > pop_sharing.csv
'''

import sys
import os
from itertools import islice
from glob import glob


# dict
site_sharing = {}

# pops
pops = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB.1KGP", "CHN.NyuWa", "CHS.1KGP",
        "CHS.NyuWa", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU",
        "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"]
for fname in glob('/home2/niuyw/project/STR/str_stat/AF/chr*.pop.txt'):
    with open(fname, 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            site = line[0] + ':' + line[1]
            count = 0
            for i,j in enumerate(range(3, 112, 4)):
                tmp = eval(line[j])
                if isinstance(tmp, tuple):
                    AC = sum(tmp)
                else:
                    AC = tmp
                if AC:
                    count += 1
            if count == 1:
                site_sharing[site] = 1 # unique
            elif count == 28:
                site_sharing[site] = 3 # all
            else:
                site_sharing[site] = 2 # shared

# each pop
lst_pops = [[0, 0, 0, 0] for i in range(28)] # tot, unique, shared, all
for i,j in enumerate(pops):
    fname = '/home2/niuyw/project/STR/str_stat/sample_stat_byPop/site_lst.%s' %(j)
    with open(fname, 'rt') as fin:
        for line in fin:
            line = line.strip()
            lst_pops[i][0] += 1
            # sharing
            count = site_sharing[line]
            lst_pops[i][count] += 1

# output
print (','.join(["pop", "tot", "Unique", "Shared", "All"]))
for i,j in enumerate(pops):
    to_w = [str(i) for i in lst_pops[i]]
    print (j + ',' + ','.join(to_w))


