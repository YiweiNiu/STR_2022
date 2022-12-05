#!/usr/bin/env python

'''
purpose: transform GWAS Catalog genes into gmt format, for ORA

usage:
cd /home2/niuyw/RefData/GWAS_Catalog
python gwasCatalog_genes_to_gmt.py
'''

import sys
from itertools import islice

# read
trait_to_genes = {}
with open('gwas_catalog_v1.0.2-associations_e100_r2020-08-26.tsv', 'rt') as fin:
    for line in islice(fin, 1, None):
        line = line.strip('\n').split('\t')
        traits = [i.strip() for i in line[34].split(',')] # MAPPED_TRAIT
        reported_genes = set([i.strip() for i in line[13].split(',')])
        mapped_genes = set([i.strip() for i in line[14].replace('-', ',').split(',')])
        genes = reported_genes | mapped_genes
        for trait in traits:
            if trait not in trait_to_genes:
                trait_to_genes[trait] = genes
            else:
                trait_to_genes[trait] = genes | trait_to_genes[trait]

# output
fout = open("gwasCatalog_genes.gmt", 'wt')
for t in trait_to_genes:
    genes = list(trait_to_genes[t])
    if len(genes) < 10:
        continue
    fout.write(t + '\tNA\t' + '\t'.join(genes) + '\n')
fout.close()



