#!/usr/bin/env python

'''
purpose:
make a long table for ExpansionHunter and str-analysis results
as the one from gnomAD: gnomAD_STR_genotypes__2022_01_20.tsv.gz

usage:
cd /home2/niuyw/project/STR/ExpansionHunter
python make_gnomAD_genotype.py
'''

import sys
from itertools import islice


VERSION = '20220413' # first run

# out
fout = open('sam6487_EH60_genotypes_20220413.tsv', 'wt')
fout.write('\t'.join(['SampleOld', 'SampleNew', 'Gender', 'Dataset', 'SuperPop', 'Population',
                      'VariantID', 'LocusId', 'RepeatType', 'ReferenceRegion', 'Motif', "MotifType",
                      'Genotype', 'GenotypeConfidenceInterval', 'ReadSupportAllele1', 'ReadSupportAllele2']) + '\n')


# read sam info
sam_info = {}
with open('/home2/niuyw/project/STR/samples/220116_sample_info.GangSTR.txt', 'rt') as fin:
    for line in islice(fin, 1, None):
        line = line.strip().split('\t')
        sex = line[3]
        if sex == "M":
            sex = 'XY'
        else:
            sex = 'XX'
        # sam_old: sam_old, sam_new, dataset, gender, super_pop, pop
        sam_info[line[1]] = [line[1], line[2], line[0], sex, line[4], line[5]]


# read repeat type
locus_2_type = {}
with open('/home2/niuyw/RefData/STR/ExpansionHunterVariantCatalog/LocusID_2_RepeatType.txt', 'rt') as fin:
    for line in islice(fin, 1, None):
        line = line.strip().split()
        locus_2_type[line[0]] = line[1]


# read ExpansionHunter results
sam_geno = {}
with open('sam6487_EH.6487_json_files.variants.tsv', 'rt') as fin:
    for line in islice(fin, 1, None):
        line = line.strip('\n').split('\t')
        sam, LocusId, VariantId = line[0], line[4], line[20]
        ReapeatType = locus_2_type[LocusId]
        # read these 9 site from str-analysis later
        if ReapeatType == 'replaced/nested':
            continue
        Motif, ReferenceRegion = line[22], line[24]
        Genotype, GenotypeConfidenceInterval = line[36], line[37]
        if '/' in Genotype:
            MotifType = 'Pathogenic/Pathogenic'
        else:
            MotifType = 'Pathogenic'
        reads_a1, reads_a2 = line[47], line[58]
        if reads_a2 == '':
            reads_a2 = 'NA'
        to_w = [VariantId, LocusId, ReapeatType, Motif, MotifType, ReferenceRegion, Genotype, GenotypeConfidenceInterval, reads_a1, reads_a2]
        fout.write('\t'.join(sam_info[sam] + to_w) + '\n')


# read str-analysis results
with open('sam6487_str_analysis.188687_json_files.tsv', 'rt') as fin:
    for line in islice(fin, 1, None):
        line = line.strip('\n').split('\t')
        sam, Filename, call = line[0], line[1], line[2]
        if call == '':
            continue
        if 'NO CALL' in call:
            continue
        MotifType = call.replace('BENIGN MOTIF', 'Benign').replace('PATHOGENIC MOTIF', 'Pathogenic').replace('MOTIF OF UNCERTAIN SIGNIFICANCE', 'Uncertain').replace(' ', '')
        LocusId = VariantId = line[29]
        ReapeatType = locus_2_type[LocusId]
        Motif, ReferenceRegion = line[6].replace(' ', ''), line[28]
        Genotype, GenotypeConfidenceInterval = line[5], line[3]
        reads_a1, reads_a2 = line[33], line[58]
        if reads_a2 == '':
            reads_a2 = 'NA'
        to_w = [VariantId, LocusId, ReapeatType, Motif, MotifType, ReferenceRegion, Genotype, GenotypeConfidenceInterval, reads_a1, reads_a2]
        fout.write('\t'.join(sam_info[sam] + to_w) + '\n')

# close
fout.close()

