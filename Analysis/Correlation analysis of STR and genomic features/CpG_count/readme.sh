#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/genomic_features_1Mb/CpG_count

# count number per 100-kb window
bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed ~/RefData/cpgIsland/hg38.220416.cpgIsland.bed > CpG.count

