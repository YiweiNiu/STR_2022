#!/usr/bin/env python

'''
purpose:
get sample pSTR stat

usage:
python3 sample_stat.py -i <xxx.vcf.gz>/<*.vcf> -o <sample_stat.csv>

sample_stat.csv
sample_old tot het period_2 period_3 period_4 period_5 period_6
sam1 xxxx
sam2 xxxx
'''

import sys
import os
import argparse
import logging
from itertools import islice
import vcfpy


VERSION = '220119' # first run

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- Get sample stat of pSTR."
    epilog = "For command line options, type: %(prog)s -h"
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('-i', '--input', nargs='+', required=True, type=str, help='Path to merged GangSTR VCF file.', metavar='', dest="input")
    argparser.add_argument('-o', default='sample_stat.csv', type=str, help='Output file.', metavar='', dest="output")
    return argparser


def run():
    """The main function
    """
    argparser = prepare_argparser()
    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    # output
    sam_stat = {} # sample stat
    # read vcfs
    for vcf in args.input:
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
        # init
        for sam in reader.header.samples.names:
            if sam not in sam_stat:
                sam_stat[sam] = [0, 0, 0, 0, 0, 0, 0] # tot,het,period_2,period_3,period_4,period_5,period_6
        # parse
        for record in reader:
            chrom, pos = record.CHROM, record.POS
            period = record.INFO['PERIOD']
            # only 2-6
            if period > 6:
                continue
            # no sex chrom
            if chrom in ["chrX", "chrY"]:
                continue
            for call in record.calls:
                sam = call.sample
                # GT:. or ./.
                if not call.called:
                    continue
                # 0/0?
                if call.is_variant:
                    sam_stat[sam][0] += 1
                else:
                    continue
                # hom/het?
                if call.is_het:
                    sam_stat[sam][1] += 1
                # period
                sam_stat[sam][period] += 1

    # output
    fout = open(args.output, 'wt')
    fout.write("sample_old,tot,het,period_2,period_3,period_4,period_5,period_6\n")
    for i in sam_stat:
        to_w = [str(i) for i in sam_stat[i]]
        fout.write(i + ',' + ','.join(to_w) + '\n')
    fout.close()


if __name__ == '__main__':
    try:
        run()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)

