#!/usr/bin/env python

'''
purpose: get distribution of alleles from NyuWa

usage:
cd /home2/niuyw/project/STR/ExpansionHunter
python get_allele_distribution.py
'''

from itertools import islice
from collections import Counter

# out dict
geno_d = {}

# fin
with open('sam6487_EH60_genotypes_20220413.tsv', 'rt') as fin:
    for line in islice(fin, 1, None):
        line = line.strip().split('\t')
        site, geno = line[6], line[12].split('/')
        reads_a1, reads_a2 = float(line[14]), line[15]
        if reads_a1 < 5:
                continue
        elif reads_a2 != 'NA' and float(reads_a2) < 5:
                continue
        if site not in geno_d:
            geno_d[site] = Counter()
        for i in geno:
            geno_d[site][i] = geno_d[site][i] + 1

# out
fout = open('sam6487.allele_distribution.txt', 'wt')
fout.write('Source\tLocus\tRepeat\tCount\n')
soure = 'NyuWa'
for site in geno_d:
    rc = geno_d[site]
    for r in rc:
        fout.write('\t'.join([soure, site, r, str(rc[r])]) + '\n')

fout.close()



