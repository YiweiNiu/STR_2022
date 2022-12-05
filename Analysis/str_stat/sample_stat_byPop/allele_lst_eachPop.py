#!/usr/bin/env python

'''
purpose: get alleles of each super pop and dataset

usage:
cd /home2/niuyw/project/STR/str_stat/sample_stat_byPop
python allele_lst_eachPop.py
'''

# dataset
datasets = {"nyuwa": "NyuWa", "1kgp": "1KGP"}
data_dir = '/Parastor300s_G30S/shiyr/gangSTR/process-all/10-loci_info'
for d in datasets:
    with open('allele_lst.%s' %(datasets[d]), 'wt') as fout:
        with open('%s/loci_info_%s.txt' %(data_dir, d), 'rt') as fin:
            for line in fin:
                line = line.strip()
                if line.startswith('chrom'):
                    continue
                line = line.split('\t')
                site = line[0] + ':' + line[1]
                acount = [i.split(':') for i in line[5].split(',')]
                for item in acount:
                    if int(item[1]) > 0:
                        allele = site + ':' + item[0]
                        fout.write(allele + '\n')

# super pop
super_pops = ["AFR", "AMR", "EAS", "EUR", "SAS"]
data_dir = '/Parastor300s_G30S/shiyr/gangSTR/process-all/14-pop'
for s in super_pops:
    with open('allele_lst.%s' %(s), 'wt') as fout:
        with open('%s/%s/%s.tab' %(data_dir, s, s), 'rt') as fin:
            for line in fin:
                line = line.strip()
                if line.startswith('chrom'):
                    continue
                line = line.split('\t')
                site = line[0] + ':' + line[1]
                acount = [i.split(':') for i in line[5].split(',')]
                for item in acount:
                    if int(item[1]) > 0:
                        allele = site + ':' + item[0]
                        fout.write(allele + '\n')


