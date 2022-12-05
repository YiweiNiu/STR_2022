#!/usr/bin/env python


'''
purpose: get a union set of variant catalog from multiple sources

usage:
cd /home2/niuyw/RefData/STR/ExpansionHunterVariantCatalog
python union_variant_catalog.py
'''

import json


# gnomAD
with open('gnomAD.variant_catalog.json', 'rt') as fin:
    gnomAD = json.load(fin)
    gnomAD_LocusId = [i['LocusId'] for i in gnomAD]
    print(len(gnomAD))

# ExpansionHunter
with open('ExpansionHunter.31.json', 'rt') as fin:
    ExpansionHunter = json.load(fin)
    for i in ExpansionHunter:
        if i['LocusId'] not in gnomAD_LocusId:
            gnomAD.append(i)
            gnomAD_LocusId.append(i['LocusId'])

# Stranger
with open('Stranger.variant_catalog.json', 'rt') as fin:
    Stranger = json.load(fin)
    for i in Stranger:
        if i['LocusId'] not in gnomAD_LocusId:
            gnomAD.append(i)
            gnomAD_LocusId.append(i['LocusId'])

# STRipy
with open('STRipy.variant_catalog.json', 'rt') as fin:
    STRipy = json.load(fin)
    for i in STRipy:
        if i['LocusId'] not in gnomAD_LocusId:
            gnomAD.append(i)
            gnomAD_LocusId.append(i['LocusId'])

print(len(gnomAD))
with open('union.variant_catalog.json', 'w') as fout:
    json.dump(gnomAD, fout)






