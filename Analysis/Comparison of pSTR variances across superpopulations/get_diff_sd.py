#!/usr/bin/env python

'''
purpose: get length var ratio (SD diff) of each STR locus for two populations

Usage:
pip install vcfpy
python3 get_diff_sd.py -a pop1.lst -b pop2.lst -v input.vcf -o pop1_vs_pop2.diff_var.txt
'''

import sys
import os
import argparse
import logging
from itertools import islice, chain
from random import shuffle
from statistics import pstdev
import vcfpy
from multiprocessing import Pool


#VERSION = '221125' # first run
VERSION = '221126' # add multiprocessing support
                   # bug fixes


# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- Get length sd diff (s1 -s2) of each STR locus between two populations."
    epilog = "For command line options, type: %(prog)s -h"
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('-t', required=True, type=int, help='Threads number.', metavar='', dest="cpu")
    argparser.add_argument('-a', required=True, type=str, help='Samples of population 1 (one sample per line, no header).', metavar='', dest="pop1")
    argparser.add_argument('-b', required=True, type=str, help='Samples of population 2 (one sample per line, no header).', metavar='', dest="pop2")
    argparser.add_argument('-v', '--vcf', required=True, type=str, help='Path to STR VCF file.', metavar='', dest="vcf")
    argparser.add_argument('-o', '--output', required=True, type=str, help='Path to output file.', metavar='', dest="output")
    return argparser


def _helper_permutation(len_all=None, n_pop1=None):
    shuffle(len_all)
    fake_pop1 = list(chain(*len_all[:n_pop1]))
    fake_pop2 = list(chain(*len_all[n_pop1:]))
    fake_sd_diff = pstdev(fake_pop1) - pstdev(fake_pop2)
    return fake_sd_diff


def get_sd_diff_p(calls=None, sam_pop1=None, sam_pop2=None, cpu=None):
    '''get SD diff and permutation p
    '''
    # repeat length
    len_all = []
    len_pop1 = []
    len_pop2 = []
    # calls
    for call in calls:
        gt = call.data.get('GT') or './.'
        if gt == './.':
            continue
        sam = call.sample
        # get length of 2 alleles
        repcn = call.data['REPCN']
        if sam in sam_pop1:
            len_pop1.extend(repcn)
            len_all.append(repcn)
        elif sam in sam_pop2:
            len_pop2.extend(repcn)
            len_all.append(repcn)
        else:
            continue # sample not to include
    # no data point
    if len(len_pop1)==0 or len(len_pop2)==0:
        return '1.0', 'NA', 'NA', 'NA'
    # sd
    sd_pop1 = pstdev(len_pop1)
    sd_pop2 = pstdev(len_pop2)
    # SD diff
    sd_diff = sd_pop1 - sd_pop2
    # permutation
    n_pop1 = int(len(len_pop1)/2)
    permutation_times = 1000
    pool = Pool(processes=cpu)
    result=[]
    for i in range(permutation_times):
        result.append(pool.apply_async(_helper_permutation, args=(len_all,n_pop1,)))
    pool.close()
    pool.join()
    # get
    fake_sd_diff_lst = []
    for i in result:
        fake_sd_diff_lst.append(i.get())
    # cmp
    if sd_diff >= 0:
        p = sum([i >= sd_diff for i in fake_sd_diff_lst])/permutation_times
    else:
        p = sum([i <= sd_diff for i in fake_sd_diff_lst])/permutation_times
    # return
    return str(p), str(sd_diff), str(sd_pop1), str(sd_pop2)


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
    cpu = args.cpu
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
    fout.write('site\tsd_pop1\tsd_pop2\tsd_diff\tpermutation.p\n')

    # get STR het
    str_het01 = {}
    with open("/home2/niuyw/project/STR/str_info/str_info.txt", 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split('\t')
            site, het = line[0], float(line[16])
            if het > 0.1:
                str_het01[site] = 0

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
        # period
        period = record.INFO['PERIOD']
        # only 2-6
        if period > 6:
            continue
        site = ":".join([chrom, str(pos)]) # chr17:32123599
        # het > 0.1
        if site not in str_het01:
            continue
        # get site-level-stat
        p, sd_diff, sd_pop1, sd_pop2 = get_sd_diff_p(record.calls, sam_pop1, sam_pop2, cpu)
        to_w = [site, sd_pop1, sd_pop2, sd_diff, p]
        # write
        fout.write('\t'.join(to_w) + '\n')
    logger.info("Complete!")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)


