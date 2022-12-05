#!/usr/bin/env python

'''
awk '/^chr[0-9]+\t/ && $4<=6' /Parastor300s_G30S/shiyr/gangSTR/1kgp/download/hg38_ver13.bed | wc -l
total sites (chr1-chr22): 765227

purpose:
get sample callrate before filtering site

usage:
python3 sample_callrate.py -i <xxx.vcf.gz>/<*.vcf> -o <sample_callrate.txt>

sample_callrate.txt:
sam1 xxxx
sam2 xxxx
'''

import sys
import os
import argparse
import logging
import vcfpy


'''
220113: first run
220114: filter site by the parameters used by dumpSTR
220119: only count STRs with period 2-6
'''

VERSION = '220119'

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- STR saturation stats."
    epilog = "For command line options, type: %(prog)s -h"
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('-i', '--input', nargs='+', required=True, type=str, help='Path to merged GangSTR VCF file.', metavar='', dest="input")
    argparser.add_argument('-o', default='sample_callrate.txt', type=str, help='Output file.', metavar='', dest="output")
    return argparser


def run():
    """The main function
    """
    argparser = prepare_argparser()
    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    # input
    vcf_lst = args.input
    samples = {}
    reasons = {} # the reasons
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
        # init
        for sam in reader.header.samples.names:
            if sam not in samples:
                samples[sam] = 0
            if sam not in reasons:
                reasons[sam] = {"Missing": 0,
                                "DPlow": 0,
                                "DPhigh": 0,
                                "SpanBoundOnly": 0,
                                "BadCI": 0}
        # parse
        for record in reader:
            chrom, pos = record.CHROM, record.POS
            # period
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
                    reasons[sam]["Missing"] += 1
                    continue
                # --gangstr-min-call-DP 20 --gangstr-max-call-DP 1000
                DP = call.data['DP']
                if DP < 20:
                    reasons[sam]["DPlow"] += 1
                    continue
                if DP > 1000:
                    reasons[sam]["DPhigh"] += 1
                    continue
                # --gangstr-filter-spanbound-only
                rcvals = eval(call.data['RC'])
                span_bound = rcvals[1] + rcvals[3]
                if span_bound == DP:
                    reasons[sam]["SpanBoundOnly"] += 1
                    continue
                # --gangstr-filter-badCI
                ml = call.data['REPCN']
                ci = call.data['REPCI']
                ci = [eval(i.replace('-', ',')) for i in ci.split(',')]
                ci1 = [i[0] for i in ci]
                ci2 = [i[1] for i in ci]
                if sum([ml[i]<j for i,j in enumerate(ci1)]) or sum([j<ml[i] for i,j in enumerate(ci2)]):
                    reasons[sam]["BadCI"] += 1
                    continue
                samples[sam] += 1
                #print (site)

    fout = open(args.output, 'wt')
    fout.write("sample_old\tcalled\tmissing\tlowDP\thighDP\tonlySpanBound\tbadCI\n")
    for sam in samples:
        to_w = [sam, str(samples[sam])] + [str(reasons[sam][i]) for i in ["Missing", "DPlow", "DPhigh", "SpanBoundOnly", "BadCI"]]
        fout.write('\t'.join(to_w) + '\n')


if __name__ == '__main__':
    try:
        run()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)

