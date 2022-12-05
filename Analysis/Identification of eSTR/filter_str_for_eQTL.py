#!/usr/bin/env python

'''
purpose: filter STRs for eQTL analysis

keep STRs:
1. loci within 500 kb of a gene expressed in our LCL data set
2. loci where at least 50 of the 445 samples had a genotype call
3. loci with heterozygosity >= 0.1
4. 

usage:
cd /home2/niuyw/project/STR/eSTR
python filter_str_for_eQTL.py --het 0.1 --numcalled 50 --tab <statSTR.out> --pair <bedtools window output> --vcf <STR.vcf>
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
    description = "%(prog)s -- Filter STRs for eQTL analysis."
    epilog = "For command line options, type: %(prog)s -h"

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('--het', default=0.1, type=float, help='Minimum heterozygosity for loci.', metavar='', dest="het")
    argparser.add_argument('--numcalled', default=50, type=int, help='Minimum number of sample called for loci.', metavar='', dest="numcalled")
    argparser.add_argument('--tab', required=True, type=str, help='Path to statSTR output.', metavar='', dest="tab")
    argparser.add_argument('--pair', required=True, type=str, help='Path to bedtools window output.', metavar='', dest="pair")
    argparser.add_argument('--vcf', required=True, type=str, help='Path to STR VCF file.', metavar='', dest="vcf")
    argparser.add_argument('-o', '--output', required=True, type=str, help='Output file prefix.', metavar='', dest="output")

    return argparser


def str_dosage(repcn_lst=None, ref=None, site=None):
    '''(Gymrek et al., 2016)
    Dosage was defined as the sum of the deviations in the STR allele lengths from the hg19 reference sequence.
    '''
    dosage_lst = []
    for x in repcn_lst:
        if not x:
            dosage_lst.append('NA')
            continue
        i,j = x
        if i:
            l = (i-ref)+(j-ref)
            dosage_lst.append(str(l))
        else:
            dosage_lst.append('NA')
    return dosage_lst


def keep_it(calls=None, ref=None, site=None):
    '''whether keep the loci?
    (Fotsing et al., 2019)
    we removed any genotypes seen fewer than three times.
    If, after filtering samples, there were fewer than three unique genotypes, the STR was excluded from analysis.
    '''
    # count
    c = Counter()
    for call in calls:
        gt = call.data.get('GT') or './.'
        c[gt] = c[gt] + 1
    # we removed any genotypes seen fewer than three times.
    gt_lst = []
    repcn_lst = []
    for call in calls:
        gt = call.data.get('GT') or './.'
        repcn = call.data['REPCN']
        if c[gt] < 3:
            gt = './.'
            repcn = [None, None]
        gt_lst.append(gt)
        repcn_lst.append(repcn)
    # count again
    c = Counter()
    for gt in gt_lst:
        if gt != "./.":
            c[gt] = c[gt] + 1
    # If, after filtering  samples, there were fewer than three unique genotypes, the STR was excluded from analysis.
    if len(c) < 3:
        return None
    else:
        dosage_lst = str_dosage(repcn_lst, ref, site)
        return dosage_lst


def main():
    """The main function
    """
    argparser = prepare_argparser()
    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    het = args.het
    numcalled = args.numcalled
    tab = args.tab
    pair = args.pair
    vcf = args.vcf
    # str to exclusion
    str_exclusion = {}
    try:
        logger.info("Read %s" %(tab))
        fin = open(tab, 'rt')
    except:
        logger.error("Read input file %s error." %(tab))
        sys.exit(1)
    else:
        for line in islice(fin, 1, None):
            line = line.strip().split('\t')
            site = line[0] + ":" + line[1]
            # het (Fotsing et al., 2019): 0.1, (Gymrek et al., 2016): 0.3
            if float(line[7]) < het:
                str_exclusion[site] = 0
                continue
            # numcalled (Fotsing et al., 2019): call rate>80%, (Gymrek et al., 2016): 50/311
            if int(line[12]) < numcalled:
                str_exclusion[site] = 0
                continue
        fin.close()
    # str-gene pair
    str_gene_pair = {}
    try:
        logger.info("Read %s" %(pair))
        fin = open(pair, 'rt')
    except:
        logger.error("Read input file %s error." %(pair))
        sys.exit(1)
    else:
        for line in fin:
            line = line.strip().split('\t')
            site = line[6] + ":" + line[7]
            # site to exclude
            if site in str_exclusion:
                continue
            # distance
            d = str(abs((int(line[2])-int(line[1])) - (int(line[8])-int(line[7]))))
            gene = line[3]
            if site not in str_gene_pair:
                str_gene_pair[site] = [(gene, d)]
            else:
                str_gene_pair[site].append((gene, d))
        fin.close()

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
    # str geno data
    fout1 = open(output + '.str_geno.csv', 'wt')
    # str-gene pair
    fout2 = open(output + '.str_gene_pair.csv', 'wt')
    fout2.write('site,gene_id,distance\n')
    # read
    try:
        logger.info("Read %s" %(vcf))
        reader = vcfpy.Reader.from_path(vcf)
    except:
        logger.error('Error when opening %s. Check it please.' %(vcf))
        sys.exit(1)
    # new str-gene pairs, since we woule continue to filter STRs here
    str_gene_pair2 = {}
    # csv, header
    fout1.write('site' + ',' + ','.join(reader.header.samples.names) + '\n')
    # parse
    for record in reader:
        chrom, pos = record.CHROM, record.POS
        site = chrom + ":" + str(pos)
        # filter
        if site not in str_gene_pair:
            continue
        # period
        period = record.INFO['PERIOD']
        # only 2-6
        if period > 6:
            continue
        # REF (length)
        REF = int(record.INFO['REF'])
        # keep it?
        dosage_lst = keep_it(record.calls, REF, site)
        if dosage_lst:
            fout1.write(site + ',' + ','.join(dosage_lst) + '\n')
            str_gene_pair2[site] = str_gene_pair[site]
    fout1.close()
    # stat
    num_str = len(str_gene_pair2)
    genes = {}
    num_pairs = 0
    for i in str_gene_pair2:
        num_pairs += len(str_gene_pair2[i])
        for g,d in str_gene_pair2[i]:
            genes[g] = 0
            fout2.write(i + ',' + g + ',' + d + '\n')
    fout2.close()
    num_genes = len(genes)
    logger.info("%s STRs." %(num_str))
    logger.info("%s genes." %(num_genes))
    logger.info("%s STR-gene pairs." %(num_pairs))
    logger.info("Complete!")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)


