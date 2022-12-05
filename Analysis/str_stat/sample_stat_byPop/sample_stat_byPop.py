#!/usr/bin/env python

'''
purpose:
get pSTR stat of pop/dataset

usage:
cd /home2/niuyw/project/STR/str_stat/sample_stat_byPop
python3 sample_stat_byPop.py > sample_stat_byPop.csv

sample_stat.csv
sample_old tot period_2 period_3 period_4 period_5 period_6
sam1 xxxx
sam2 xxxx
'''

import sys
import os
from itertools import islice

VERSION = '220121' # first run

# read site period
site_period = {}
with open('/home2/niuyw/project/STR/str_stat/str_info/str_info.csv', 'rt') as fin:
    for line in islice(fin, 1, None):
        line = line.strip().split(',')
        site = line[0] + ':' + line[1]
        period = int(line[4])
        site_period[site] = period

# datasets
datasets = ["Merged", "1KGP", "NyuWa"]
lst_dataset = [[i, 0, 0, 0, 0, 0, 0] for i in datasets] # dataset, tot, period_2, period_3, period_4, period_5, period_6

for i,j in enumerate(datasets):
    fname = 'site_lst.%s' %(j)
    with open(fname, 'rt') as fin:
        for line in fin:
            line = line.strip()
            period = site_period[line]
            lst_dataset[i][1] += 1
            lst_dataset[i][period] += 1

# super pop
super_pops = ["AFR", "AMR", "EAS", "EUR", "SAS"]
lst_spop = [[i, 0, 0, 0, 0, 0, 0] for i in super_pops]
for i,j in enumerate(super_pops):
    fname = 'site_lst.%s' %(j)
    with open(fname, 'rt') as fin:
        for line in fin:
            line = line.strip()
            period = site_period[line]
            lst_spop[i][1] += 1
            lst_spop[i][period] += 1

# pop
pops = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB.1KGP", "CHN.NyuWa", "CHS.1KGP",
        "CHS.NyuWa", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU",
        "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"]
lst_pop = [[i, 0, 0, 0, 0, 0, 0] for i in pops]
for i,j in enumerate(pops):
    fname = 'site_lst.%s' %(j)
    with open(fname, 'rt') as fin:
        for line in fin:
            line = line.strip()
            period = site_period[line]
            lst_pop[i][1] += 1
            lst_pop[i][period] += 1

# print
print (','.join(["dataset", "tot", "period_2", "period_3", "period_4", "period_5", "period_6"]))
for i in lst_dataset:
    to_w = [str(x) for x in i]
    print (','.join(to_w))
for i in lst_spop:
    to_w = [str(x) for x in i]
    print (','.join(to_w))
for i in lst_pop:
    to_w = [str(x) for x in i]
    print (','.join(to_w))

