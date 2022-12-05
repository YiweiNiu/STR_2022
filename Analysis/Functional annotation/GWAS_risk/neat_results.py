#!/usr/bin/env python


'''
purpose: neat the results from LD of GWAS SNP and STRs

usage:
cd /home2/niuyw/project/STR/GWAS_risk
python neat_results.py
'''

import sys
from glob import glob


def shear(traits=None, pvalues=None):
    t_p_dict = {}
    for idx,t in enumerate(traits):
        if t not in t_p_dict:
            t_p_dict[t] = float(pvalues[idx])
        else:
            p = float(pvalues[idx])
            if p < t_p_dict[t]:
                t_p_dict[t] = float(p)
    t_lst, p_lst = [], []
    for t in t_p_dict:
        t_lst.append(t)
        p_lst.append(str(t_p_dict[t]))
    return t_lst, p_lst


def read_snp_info(fname=None):
    snp_info = {}
    with open(fname, 'rt') as fin:
        for line in fin:
            line = line.strip().split('\t')
            snp = line[0] + ":" + line[1]
            rs_id = line[2].split('-')[0]
            traits, pvalues = shear(line[4].split('@'), line[5].split('@'))
            snp_info[snp] = [rs_id, traits, pvalues]
    return snp_info


def read_str_snp_ld():
    str_snp_ld = {}
    for fname in glob('LD/*.ld.txt'):
        with open(fname, 'rt') as fin:
            for line in fin:
                line = line.strip().split()
                if line[6] == 'nan':
                    continue
                ld = float(line[6])
                if ld < 0.7:
                    continue
                STR = line[3] + ':' + line[4]
                snp = line[0] + ':' + line[2]
                if STR not in str_snp_ld:
                    str_snp_ld[STR] = [(snp, str(ld))]
                else:
                    str_snp_ld[STR].append((snp, str(ld)))
    return str_snp_ld


def main():
    snp_info_fname = 'gwas_snps.txt'
    snp_info = read_snp_info(snp_info_fname)
    str_snp_ld = read_str_snp_ld()
    fout = open('STR_anno_byGWASCatalog.txt', 'wt')
    fout.write('\t'.join(["STR", "SNP_POS", "SNP_ID", "Trait", "PValue", "LD"]) + '\n')
    for s in str_snp_ld:
        for item in str_snp_ld[s]:
            snp, ld = item
            rs_id, traits, pvalues = snp_info[snp]
            for idx,t in enumerate(traits):
                fout.write('\t'.join([s, snp, rs_id, t, pvalues[idx], ld]) + '\n')
    fout.close()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted me! ;-) Bye!")
        sys.exit(1)

