#!/usr/bin/env python

'''
purpose: collect the results of sample_callrate.py

usage:
cd /home2/niuyw/project/STR/str_stat/sample_stat
python sample_stat_collect.py chr*.csv > sample_stat.csv
'''

import sys
from itertools import islice

sam_dict = {}
for fname in sys.argv[1:]:
    with open(fname, 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split(',')
            sam = line[0]
            called = [int(i) for i in line[1:]]
            if sam not in sam_dict:
                sam_dict[sam] = called
            else:
                sam_dict[sam] = [sum(x) for x in zip(sam_dict[sam], called)]

# write
print("sample_old,tot,het,period_2,period_3,period_4,period_5,period_6")
for sam in sam_dict:
    print (sam + ',' + ','.join([str(i) for i in sam_dict[sam]]))


