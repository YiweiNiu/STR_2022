#!/usr/bin/env python

'''
purpose: calculate STR-SNP LDs

usage:
cd /home2/niuyw/project/STR/GWAS_risk
python get_SNP_STR_LD.py
'''

import sys
import os
import argparse
import logging
import pandas as pd
from multiprocessing import Pool
from functools import partial


VERSION = '220323' # first run

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
    argparser.add_argument('--cpu', default=4, type=int, help='CPU numbers.', metavar='', dest="cpu")
    argparser.add_argument('--pair', required=True, type=str, help='Path to bedtools window output.', metavar='', dest="pair")
    argparser.add_argument('--chrom', required=True, type=str, help='Chrom to handle.', metavar='', dest="chrom")
    argparser.add_argument('--snp', required=True, type=str, help='Path to SNP dosage matrix.', metavar='', dest="snp")
    argparser.add_argument('--str', required=True, type=str, help='Path to STR dosage matrix.', metavar='', dest="str")
    argparser.add_argument('--output', required=True, type=str, help='Output file prefix.', metavar='', dest="output")
    return argparser


def read_snp_tr_pair(fname=None, chrom=None):
    snp_str_pair = {}
    STRs = {}
    with open(fname, 'rt') as fin:
        for line in fin:
            line = line.strip('\n').split('\t')
            if line[0] != chrom:
                continue
            snp = line[0] + ":" + line[2]
            tr = line[3] + ":" + line[4]
            if snp not in snp_str_pair:
                snp_str_pair[snp] = [tr]
            else:
                snp_str_pair[snp].append(tr)
            STRs[tr] = 0
    return snp_str_pair, STRs


def cal_snp_str_ld_helper(x=None, y=None):
    if isinstance(y, pd.DataFrame):
        cor_lst = []
        for ix,col in y.iteritems():
            cor = x.corr(col)
            cor_lst.append(cor)
        cor = max(cor_lst)
        ld = cor*cor
    elif isinstance(y, pd.Series):
        cor = x.corr(y)
        ld = cor*cor
    else:
        ld = 'nan'
    return str(ld)


def cal_snp_str_ld(df_snp = None, df_str = None, snp=None, tr=None):
    # str
    x = df_str[tr]
    # snp
    snp_chrom, snp_pos = snp.split(':')
    snp2 = snp_chrom + ":" + str(int(snp_pos)-1) # also possible
    if snp in df_snp.columns:
        y = df_snp[snp]
        ld = cal_snp_str_ld_helper(x, y)
    elif snp2 in df_snp.columns:
        y = df_snp[snp2]
        ld = cal_snp_str_ld_helper(x, y)
    else:
        ld = 'nan'
    return [snp, tr, ld]


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
    # read SNP-STR pairs
    snp_str_pair, STRs = read_snp_tr_pair(args.pair, chrom)
    logger.info("Load SNP-STR pairs done.")
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
    # get LD
    fout = open('%s.ld.txt' %(output), 'wt')
    logger.info("Computing LDs...")
    results = []
    # pool
    p = Pool(args.cpu)
    for snp in snp_str_pair:
        for tr in snp_str_pair[snp]:
            res = p.apply_async(cal_snp_str_ld, args=(df_snp, df_str, snp, tr,))
            results.append(res)
    p.close()
    p.join()
    # get results
    for res in results:
        snp, tr, ld = res.get()
        snp_chrom, snp_pos = snp.split(':')
        tr_chrom, tr_pos = tr.split(':')
        fout.write('\t'.join([snp_chrom, snp_pos, tr_chrom, tr_pos, ld]) + '\n')
    fout.close()
    logger.info("All done!")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)


