#!/usr/bin/env python

'''
purpose:
get pSTR lst of each pop

usage:
cd /home2/niuyw/project/STR/str_stat/sample_stat_byPop
python3 site_lst_eachPop.py
'''

import sys
import os
from itertools import islice
from glob import glob

# dataset
datasets = ["1KGP", "NyuWa"]
lst_dataset = [{}, {}] # "1KGP", "NyuWa"
for fname in glob('/home2/niuyw/project/STR/str_stat/AF/chr*.dataset.txt'):
    with open(fname, 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            site = line[0] + ':' + line[1]
            for i,j in enumerate([7, 11]):
                tmp = eval(line[j])
                if isinstance(tmp, tuple):
                    AC = sum(tmp)
                else:
                    AC = tmp
                if AC:
                    lst_dataset[i][site] = 0
for i,j in enumerate(datasets):
    with open('site_lst.%s' %(j), 'wt') as fout:
        for s in lst_dataset[i]:
            fout.write(s + '\n')

# super pop
super_pops = ["AFR", "AMR", "EAS", "EUR", "SAS"]
lst_spop = [{}, {}, {}, {}, {}, {}] # AFR, AMR, EAS, EUR, SAS
for fname in glob('/home2/niuyw/project/STR/str_stat/AF/chr*.superPop.txt'):
    with open(fname, 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            site = line[0] + ':' + line[1]
            for i,j in enumerate([3, 7, 11, 15, 19]):
                tmp = eval(line[j])
                if isinstance(tmp, tuple):
                    AC = sum(tmp)
                else:
                    AC = tmp
                if AC:
                    lst_spop[i][site] = 0
for i,j in enumerate(super_pops):
    with open('site_lst.%s' %(j), 'wt') as fout:
        for s in lst_spop[i]:
            fout.write(s + '\n')

# pops
pops = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB.1KGP", "CHN.NyuWa", "CHS.1KGP",
        "CHS.NyuWa", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU",
        "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"]
lst_pops = [dict() for i in range(28)]
for fname in glob('/home2/niuyw/project/STR/str_stat/AF/chr*.pop.txt'):
    with open(fname, 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            site = line[0] + ':' + line[1]
            for i,j in enumerate(range(3, 112, 4)):
                tmp = eval(line[j])
                if isinstance(tmp, tuple):
                    AC = sum(tmp)
                else:
                    AC = tmp
                if AC:
                    lst_pops[i][site] = 0
for i,j in enumerate(pops):
    with open('site_lst.%s' %(j), 'wt') as fout:
        for s in lst_pops[i]:
            fout.write(s + '\n')


