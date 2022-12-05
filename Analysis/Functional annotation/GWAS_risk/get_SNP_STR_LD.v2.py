#!/usr/bin/env python

'''
purpose: calculate STR-SNP LDs

usage:
cd /home2/niuyw/project/STR/GWAS_risk
chrom=chr21
python get_SNP_STR_LD.v2.py --pair gwas_snps.window250k_pSTRs.txt --chrom $chrom --snp SNPs/$chrom.raw --str STRs/$chrom.str_dosage.csv --output LD/$chrom
'''

import sys
import os
import argparse
import logging
from collections import defaultdict
from copy import deepcopy
import pandas as pd
from numpy import nan
from math import isnan


#VERSION = '220323' # first run
VERSION = '220325' # use pandas to calculate correlation

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- Calculate LD between SNPs and STRs."
    epilog = "For command line options, type: %(prog)s -h"
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('--pair', required=True, type=str, help='Path to bedtools window output.', metavar='', dest="pair")
    argparser.add_argument('--chrom', required=True, type=str, help='Chrom to handle.', metavar='', dest="chrom")
    argparser.add_argument('--snp', required=True, type=str, help='Path to SNP dosage matrix.', metavar='', dest="snp")
    argparser.add_argument('--str', required=True, type=str, help='Path to STR dosage matrix.', metavar='', dest="str")
    argparser.add_argument('--output', required=True, type=str, help='Output file prefix.', metavar='', dest="output")
    return argparser


def get_snp_tr_ld(fname=None, chrom=None, dict_cor=None, renamer=None, fout=None):
    with open(fname, 'rt') as fin:
        for line in fin:
            line = line.strip('\n')
            items = line.split('\t')
            if items[0] != chrom:
                continue
            snp = items[0] + ":" + items[2]
            snp2 = items[0] + ":" + items[1]
            tr = items[3] + ":" + items[4]
            if snp in dict_cor:
                cor = dict_cor[snp][tr]
            elif snp2 in dict_cor:
                cor = dict_cor[snp2][tr]
            elif snp in renamer:
                tmp_lst = []
                for s in renamer[snp]:
                    cor = dict_cor[s][tr]
                    if not isnan(cor):
                        tmp_lst.append(cor)
                if tmp_lst:
                    cor = max(tmp_lst)
                else:
                    cor = nan
            elif snp2 in renamer:
                tmp_lst = []
                for s in renamer[snp2]:
                    cor = dict_cor[s][tr]
                    if not isnan(cor):
                        tmp_lst.append(cor)
                if tmp_lst:
                    cor = max(tmp_lst)
                else:
                    cor = nan
            else:
                cor = nan
            ld = str(cor*cor)
            fout.write(line + '\t' + ld + '\n')


def main():
    """The main function
    """
    argparser = prepare_argparser()
    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    chrom = args.chrom
    # output
    output = args.output
    output_path = os.path.abspath(os.path.dirname(output))
    if not os.path.exists(output_path): os.mkdir(output_path)
    logger.info("Start analysis.")
    # Load SNP dosage
    df_snp = df_snp = pd.read_csv(args.snp, sep=' ')
    df_snp.drop(df_snp.columns[[0, 2, 3, 4, 5]], axis=1, inplace=True)
    logger.info("Load SNP dosage done.")
    # Load STR dosage
    df_str = pd.read_csv(args.str, sep=',')
    df_str = df_str.set_index('site')
    df_str = df_str.T
    logger.info("Load STR dosage done.")
    # sort STR table by SNP table
    df_str = df_str.reindex(index=df_snp['IID'])
    df_snp = df_snp.set_index('IID')
    # change names of SNP
    x = df_snp.columns.to_list()
    x = {i:'chr' + ':'.join(i.split(':')[:2]) for i in x}
    df_snp.rename(columns=x, inplace=True)
    # duplicate columns
    renamer = defaultdict()
    for column_name in df_snp.columns[df_snp.columns.duplicated(keep=False)].tolist():
        if column_name not in renamer:
            renamer[column_name] = [column_name+'_0']
        else:
            renamer[column_name].append(column_name +'_'+str(len(renamer[column_name])))
    renamer2 = deepcopy(renamer)
    df_snp.rename(
        columns=lambda column_name: renamer[column_name].pop(0)
        if column_name in renamer 
        else column_name,
        inplace=True
    )
    logger.info("Neat data done.")
    # get LD
    fout = open('%s.ld.txt' %(output), 'wt')
    logger.info("Computing LDs...")
    df_cor = pd.concat([df_snp, df_str], axis=1, keys=['df_snp', 'df_str']).corr().loc['df_snp', 'df_str']
    dict_cor = df_cor.to_dict("index")
    get_snp_tr_ld(args.pair, chrom, dict_cor, renamer2, fout)
    fout.close()
    logger.info("All done!")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)


