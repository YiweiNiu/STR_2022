#!/usr/bin/env python

'''
purpose: get SNPs with p<=5E-8 from GWAS Catalog

usage:
cd /home2/niuyw/project/STR/GWAS_risk

'''

import sys
import os


def parse_gwas(file_dir):
    gwas_file = open(file_dir,'r')
    return_dict = {str(i):{} for i in range(1,23)}
    line = gwas_file.readline()
    line = gwas_file.readline()
    while line:
        li = line.split('\t')
        snp_id, chrom, pos, trait_reported, trait, pvalue = (li[20]), li[11], li[12], li[7], li[34], float(li[27])
        if pvalue>5e-8 or (chrom not in return_dict):
            line = gwas_file.readline()
            continue
        else:
            if pos not in return_dict[chrom]:
                return_dict[chrom][pos] = {'snp_id':snp_id, 'trait_reported':[trait_reported], 'trait':[trait], 'pvalue':[str(pvalue)]}
            else:
                return_dict[chrom][pos]['snp_id'] = snp_id
                return_dict[chrom][pos]['trait_reported'].append(trait_reported)
                return_dict[chrom][pos]['trait'].append(trait)
                return_dict[chrom][pos]['pvalue'].append(str(pvalue))
        line = gwas_file.readline()
    gwas_file.close()
    return return_dict


def main():
    gwas_snps = parse_gwas('/home2/niuyw/RefData/GWAS_Catalog/gwas_catalog_v1.0.2-associations_e100_r2020-08-26.tsv')
    for chrom in gwas_snps:
        for pos in gwas_snps[chrom]:
            tmp = gwas_snps[chrom][pos]
            trait_reported = '@'.join(tmp['trait_reported'])
            trait = '@'.join(tmp['trait'])
            pvalue = '@'.join(tmp['pvalue'])
            print ('\t'.join(['chr'+chrom, pos, tmp['snp_id'], trait_reported, trait, pvalue]))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)

