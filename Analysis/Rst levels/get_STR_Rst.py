#!/usr/bin/env python

'''
purpose: get Rst of each STR locus for two populations

Usage:
pip install vcfpy
python3 get_STR_Rst.py -a pop1.lst -b pop2.lst -v input.vcf -o pop1_vs_pop2.Rst
python get_STR_Rst.py -a lst.CHN.old -b lst.CHS.old -v /Parastor300s_G30S/shiyr/v7/v7-dump-top100-sort.vcf -o test.txt

output:
site Rst
chr17:32123599:tttc:4 0.01
'''

import sys
import os
import argparse
import logging
import vcfpy
from statistics import variance

'''
220109: test /Parastor300s_G30S/shiyr/v7/v7-dump-top100-sort.vcf
'''

#VERSION = '220110'
VERSION = '220409' # correct error when calculating Rst

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- Get Rst of each STR locus between two populations."
    epilog = "For command line options, type: %(prog)s -h"

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('-a', required=True, type=str, help='Samples of population 1 (one sample per line, no header).', metavar='', dest="pop1")
    argparser.add_argument('-b', required=True, type=str, help='Samples of population 2 (one sample per line, no header).', metavar='', dest="pop2")
    argparser.add_argument('-v', '--vcf', required=True, type=str, help='Path to STR VCF file.', metavar='', dest="vcf")
    argparser.add_argument('-o', '--output', required=True, type=str, help='Path to output file.', metavar='', dest="output")

    return argparser


def get_rst(calls=None, sam_pop1=None, sam_pop2=None):
    '''get Rst
    '''
    # repeat length
    len_pop1 = []
    len_pop2 = []

    for call in calls:
        gt = call.data.get('GT') or './.'
        if gt == './.': continue
        sam = call.sample
        # get length of 2 alleles
        if sam in sam_pop1:
            len_pop1.append(sum(call.data['REPCN']))
        elif sam in sam_pop2:
            len_pop2.append(sum(call.data['REPCN']))
        else:
            continue # sample not to include
    # no data point
    if len(len_pop1)<2 or len(len_pop2)<2:
        return '0', str(len(len_pop1)), str(len(len_pop2))
    S_w = (variance(len_pop1) + variance(len_pop2))/2.0
    S_bar = variance(len_pop1+len_pop2)
    if S_bar != 0:
        Rst = (S_bar - S_w)/S_bar
    else:
        Rst = 0
    # Right?
    # Negative Fst values should be effectively seen as zero values.
    # A zero value for Fst means that there is no genetic subdivision between the populations considered
    if Rst < 0:
        Rst = 0
    return str(Rst), str(len(len_pop1)), str(len(len_pop2))


def main():
    """The main function
    """
    argparser = prepare_argparser()
    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    pop1 = args.pop1
    pop2 = args.pop2
    vcf = args.vcf
    # sample info
    sam_pop1 = {}
    try:
        logger.info("Read %s" %(pop1))
        fin = open(pop1, 'rt')
    except:
        logger.error("Read input file %s error." %(pop1))
        sys.exit(1)
    else:
        for line in fin:
            line = line.strip()
            sam_pop1[line] = 0
        fin.close()
        logger.info("There are %s samples in %s." %(len(sam_pop1), pop1))
    sam_pop2 = {}
    try:
        logger.info("Read %s" %(pop2))
        fin = open(pop2, 'rt')
    except:
        logger.error("Read input file %s error." %(pop2))
        sys.exit(1)
    else:
        for line in fin:
            line = line.strip()
            sam_pop2[line] = 0
        fin.close()
        logger.info("There are %s samples in %s." %(len(sam_pop2), pop2))
    # check olp of sam_pop1 and sam_pop2
    tmp_len = len(set(sam_pop1) & set(sam_pop2))
    if tmp_len > 0:
        logger.error("Some %s samples are in both %s and %s." %(tmp_len, pop1, pop2))
        sys.exit(1)

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
    fout = open(output, 'wt')
    fout.write('site\talleleNum_pop1\talleleNum_po2\tRst\n')

    # read
    try:
        logger.info("Read %s" %(vcf))
        reader = vcfpy.Reader.from_path(vcf)
    except:
        logger.error('Error when opening %s. Check it please.' %(vcf))
        sys.exit(1)

    # parse
    for record in reader:
        chrom, pos = record.CHROM, record.POS
        site = chrom + ":" + str(pos) # chr22:10510344
        # get site-level-stat
        Rst, len_pop1, len_pop2 = get_rst(record.calls, sam_pop1, sam_pop2)
        to_w = [site, len_pop1, len_pop2, Rst]
        # write
        fout.write('\t'.join(to_w) + '\n')
    logger.info("Complete!")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)


