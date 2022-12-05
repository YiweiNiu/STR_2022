#!/usr/bin/env python3

'''
purpose:
neat gtcheck out (modified from tengxy)

usage: python3 norm.py -i <xxx.txt>/<*.txt> -o <str_diversity.txt>

output
population  str_diversity   pair_count
CHN 963.8702799897252   105111
CHS 953.1584376430052   185745
GBR 922.0587601078167   5565
FIN 895.2804029304029   5460
CHS.1KGP    868.3357464607465   6216
'''

import sys
import os
import argparse

'''change log
220112: test run
'''

VERSION = '220112'

def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- Neat gtcheck."
    epilog = "For command line options, type: %(prog)s -h"
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('-i', '--input', nargs='+', required=True, type=str, help='Path to outputs of bcftools gtcheck.', metavar='', dest="input")
    argparser.add_argument('-o', default='str_diversity.txt', type=str, help='Output file.', metavar='', dest="output")
    return argparser


def parse_sample_detail(sample_file_path=None):
    sample_file = open(sample_file_path, 'rt')
    line = sample_file.readline()
    line = sample_file.readline()
    race_sample = {}
    while line:
        line = line.strip('\n')
        li = line.split('\t')
        pop = li[5]
        if pop == 'NA':
            line = sample_file.readline()
            continue
        if pop not in race_sample:
            race_sample[pop] = set()
        race_sample[pop].add(li[1])
        line = sample_file.readline()
    sample_file.close()
    return race_sample


def process_file(input_files=None, race_sample=None):
    # init dict
    diversity_dict = {}
    count_dict = {}
    for race in race_sample:
        count_dict[race] = 0
        diversity_dict[race] = 0
    # read files
    for fname in input_files:
        input_file = open(fname, 'r')
        for line in input_file:
            if not line.startswith('CN'):
                continue
            line = line.strip()
            li = line.split('\t')
            for race in race_sample:
                if li[4] in race_sample[race] and li[5] in race_sample[race]:
                    diversity_dict[race] += int(li[1])
                    count_dict[race] += 1
        input_file.close()
    # update
    for race in race_sample:
        if count_dict[race] != 0:
            diversity_dict[race] = float(diversity_dict[race]) / float(count_dict[race])
        else:
            diversity_dict.pop(race)
            count_dict.pop(race)
    return diversity_dict, count_dict


if __name__ == '__main__':
    argparser = prepare_argparser()
    args = argparser.parse_args()

    race_sample = parse_sample_detail('/home2/niuyw/project/STR/samples/220116_sample_info.GangSTR.txt')

    input_files = args.input
    output_file = args.output
    average_diversity, count_dict = process_file(input_files, race_sample)

    with open(output_file, 'wt') as out_file:
        out_file.write('population\tstr_diversity\tpair_count\n')
        for race in average_diversity:
            out_file.write(race+'\t'+str(average_diversity[race])+'\t'+str(count_dict[race])+'\n')



