#!/usr/bin/env python

'''
purpose: STR to matrix

how:
LD between STRs and SNPs was computed by taking the squared Pearson correlation between STR lengths and SNP dosages in GTEx samples for each STR–SNP pair. (Fotsing et al., 2019)

Pairwise-LD between the eSTR and eSNPs was estimated using the Pearson correlation between SNP dosages (0, 1, or 2) and STR dosages (sum of the two repeat allele lengths). (Saini et al., 2018)

usage:
cd /home2/niuyw/project/STR/eSTR
python str_to_mat.py --het 0.1 --numcalled 50 --tab <statSTR.out> --pair <bedtools window output> --vcf <STR.vcf>
'''

import sys
import os
import argparse
from itertools import islice
from collections import Counter
import logging
import vcfpy


VERSION = '220309' # first run

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- STR to matrix."
    epilog = "For command line options, type: %(prog)s -h"

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('-v', required=True, type=str, help='Path to STR VCF file.', metavar='', dest="vcf")
    argparser.add_argument('-o', required=True, type=str, help='Output file prefix.', metavar='', dest="output")
    return argparser


def str_dosage_helper(repcn_lst=None):
    '''(Saini et al., 2018)
    sum of the two repeat allele lengths
    '''
    dosage_lst = []
    for x in repcn_lst:
        if not x:
            dosage_lst.append('NA')
            continue
        i,j = x
        if i:
            l = i+j
            dosage_lst.append(str(l))
        else:
            dosage_lst.append('NA')
    return dosage_lst


def str_dosage(calls=None):
    '''whether keep the loci?
    (Fotsing et al., 2019)
    STR genotypes seen less than three times were filtered from LD calculations.
    '''
    # count
    c = Counter()
    for call in calls:
        gt = call.data.get('GT') or './.'
        c[gt] = c[gt] + 1
    # we removed any genotypes seen fewer than three times.
    repcn_lst = []
    for call in calls:
        gt = call.data.get('GT') or './.'
        repcn = call.data['REPCN']
        if c[gt] < 3:
            gt = './.'
            repcn = [None, None]
        repcn_lst.append(repcn)
    dosage_lst = str_dosage_helper(repcn_lst)
    return dosage_lst


def main():
    """The main function
    """
    argparser = prepare_argparser()
    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    vcf = args.vcf
    # vcf
    if  os.path.exists(vcf) and os.path.getsize(vcf) > 0:
        pass
    else:
        logger.error("Read input file %s error." %(vcf))
        sys.exit(1)
    # output
    output = args.output
    output_path = os.path.abspath(os.path.dirname(output))
    if not os.path.exists(output_path): os.mkdir(output_path)
    # str dosage data
    fout = open(output + '.str_dosage.csv', 'wt')
    # read
    try:
        logger.info("Read %s" %(vcf))
        reader = vcfpy.Reader.from_path(vcf)
    except:
        logger.error('Error when opening %s. Check it please.' %(vcf))
        sys.exit(1)
    # csv, header
    fout.write('site' + ',' + ','.join(reader.header.samples.names) + '\n')
    # parse
    for record in reader:
        chrom, pos = record.CHROM, record.POS
        site = chrom + ":" + str(pos)
        # only 2-6
        if record.INFO['PERIOD'] > 6:
            continue
        # dosage
        dosage_lst = str_dosage(record.calls)
        fout.write(site + ',' + ','.join(dosage_lst) + '\n')
    fout.close()
    logger.info("Complete!")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)


