#!/usr/bin/env python

'''
purpose: convert vcf of STR to matrix

Usage:
pip install vcfpy
python3 str_vcf2mat.py -f 0.95 -t [indicator, length] -i input.vcf.gz -o output.mat

output:
locus sam1 sam2
chr1:wew?: 1 0
chr2:wewe:wew 0 9
'''

import sys
import os
import argparse
import logging
import vcfpy

# change log
#VERSION = '220113' # test /Parastor300s_G30S/shiyr/v7/v7-dump-top100-sort.vcf
VERSION = '220115' # filter by MAF


# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- STR vcf to matrix."
    epilog = "For command line options, type: %(prog)s -h"
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument("-t", choices=["indicator", "length"], default="indicator", type=str, help="Which type of matrix?", metavar="", dest="mat_type")
    argparser.add_argument('-f', default=0.95, type=float, help='Maximum major allele frequency for site selection.', metavar='', dest="freq")
    argparser.add_argument('-i', '--input', required=True, type=str, help='Path to merged GangSTR VCF file.', metavar='', dest="input")
    argparser.add_argument('-o', required=True, type=str, help='Output file.', metavar='', dest="output")
    return argparser


def run():
    """The main function
    """
    argparser = prepare_argparser()
    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    # matrix type
    mat_type = args.mat_type
    # maf
    freq = args.freq
    # output
    output = args.output
    output_path = os.path.abspath(os.path.dirname(output))
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    fout = open(output, 'wt')
    fout.write('Allele,')
    vcf = args.input
    # input
    if  os.path.exists(vcf) and os.path.getsize(vcf) > 0:
        pass
    else:
        "Input file error: empty file or file not exist"
    # read
    try:
        reader = vcfpy.Reader.from_path(vcf)
    except:
        print ('Error when opening %s. Check it please.' %(vcf))
        sys.exit(1)
    # samples
    samples = reader.header.samples.names
    sam_num = len(samples)
    fout.write(','.join(samples) + '\n')
    # parse
    for record in reader:
        chrom, pos = record.CHROM, record.POS
        # no alt
        if not record.ALT:
            continue
        # period
        period = record.INFO['PERIOD']
        # only 2-6
        if period > 6:
            continue
        # RU, REF (length)
        RU, REF = record.INFO['RU'], int(record.INFO['REF'])
        site = ":".join([chrom, str(pos), RU, str(REF)]) # chr22:10510344:taa:5
        # filter by maf (major allele frequency)
        maf = max([1-sum(record.INFO['AF']), max(record.INFO['AF'])])
        if maf >= freq:
            continue
        # to mat
        if mat_type == "indicator":
            # init output
            alt_num = len(record.ALT)
            geno_lst = [[0]*sam_num for i in range(alt_num)] # 0:alt num
            for call in record.calls:
                sam_idx = samples.index(call.sample)
                if call.called:
                    for i in call.gt_alleles:
                        if i == 0:
                            continue
                        else:
                            geno_lst[i-1][sam_idx] += 1
                else:
                    for lst in geno_lst:
                        lst[sam_idx] = 9 # missing
            # ouput
            for idx,lst in enumerate(geno_lst):
                name = site + ':a%s' %(idx)
                lst = [str(i) for i in lst]
                fout.write(name + ',' + ','.join(lst) + '\n')
        else:
            # init output
            geno_lst = [0]*sam_num
            for call in record.calls:
                sam_idx = samples.index(call.sample)
                if call.called:
                    REPCN = call.data['REPCN']
                    geno_lst[sam_idx] = str(sum(REPCN))
                else:
                    geno_lst[sam_idx] = 'NA'
            # ouput
            fout.write(site + ',' + ','.join(geno_lst) + '\n')
    fout.close()


if __name__ == '__main__':
    try:
        run()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)


