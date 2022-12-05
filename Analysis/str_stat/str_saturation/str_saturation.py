#!/usr/bin/env python

'''
purpose:
将 sample 顺序随机打乱 20 次，统计每次 100 个sample
每增加 100 个sample，计数各种变异的个数
major allele frequency
1. polymorphic STR -- pSTR (non ./., 0/0)
2. MAF (0, 0.99] common -- commonSTR
3. MAF (0.99, 1] rare -- rareSTR

usage: python3 str_saturation.py -t 3983 -i <xxx.vcf.gz>/<*.vcf> -o <str_saturation>

str_saturation.txt:
xxx xxx xxx
'''

import sys
import os
import argparse
import vcfpy
import random
import numpy

# change log
#VERSION = '220112' # first run
#VERSION = '220119' # only count STRs with period 2-6
VERSION = '220509' # use major allele frequency to replace minor allele frequency


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- STR saturation stats."
    epilog = "For command line options, type: %(prog)s -h"
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('-i', '--input', nargs='+', required=True, type=str, help='Path to merged GangSTR VCF file.', metavar='', dest="input")
    argparser.add_argument('-t', '--total', required=True, type=int, help='Total number of samples.', metavar='', dest="total")
    argparser.add_argument('-o', default='str_saturation', type=str, help='Prefix output file.', metavar='', dest="output")
    return argparser


def run():
    """The main function
    """
    argparser = prepare_argparser()

    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)

    args = argparser.parse_args()

    # total samples
    tot_n = args.total
    row_n = int(tot_n/100) + 1
    dic_sta={'pSTR': numpy.zeros([row_n, 20]),
             'commonSTR': numpy.zeros([row_n, 20]),
             'rareSTR': numpy.zeros([row_n, 20]),
             'veryrareSTR': numpy.zeros([row_n, 20])}

    # input
    vcf_lst = args.input
    for vcf in vcf_lst:
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

        # prepare
        sample_indexs = numpy.arange(0, tot_n)
        sample_group = [None] * tot_n
        for i in sample_indexs:
            sample_group[i] = []
        for k in range(20):
            random.shuffle(sample_indexs)
            for i, j in enumerate(sample_indexs):
                i2 = int(i/100)
                p = i2, k
                sample_group[j].append(p)

        # parse
        for record in reader:
            chrom, pos = record.CHROM, record.POS
            # period
            period = record.INFO['PERIOD']
            # only 2-6
            if period > 6:
                continue
            maf = record.INFO['MAF']
            if maf >= 0.01:
                var = 'commonSTR'
            elif maf < 0.001:
                var = 'veryrareSTR'
            else:
                var = 'rareSTR'
            tmp = numpy.zeros([row_n, 20])
            for j, call in enumerate(record.calls):
                if call.data['GT'] in ('0/0', './.'): continue
                for p in sample_group[j]:
                    i, k = p
                    if tmp[i, k] == 1: continue
                    for i2 in range(i, row_n):
                      tmp[i2, k] = 1

            dic_sta[var] += tmp
            dic_sta['pSTR'] += tmp

    output = args.output
    file_w1 = '%s.pSTR.csv' %(output)
    numpy.savetxt(file_w1, dic_sta['pSTR'], fmt="%d", delimiter=",")
    file_w2 = '%s.commonSTR.csv' %(output)
    numpy.savetxt(file_w2, dic_sta['commonSTR'], fmt="%d", delimiter=",")
    file_w3 = '%s.rareSTR.csv' %(output)
    numpy.savetxt(file_w3, dic_sta['rareSTR'], fmt="%d", delimiter=",")
    file_w4 = '%s.veryrareSTR.csv' %(output)
    numpy.savetxt(file_w4, dic_sta['veryrareSTR'], fmt="%d", delimiter=",")


if __name__ == '__main__':
    try:
        run()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)

