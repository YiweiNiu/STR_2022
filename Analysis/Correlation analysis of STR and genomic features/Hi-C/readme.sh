#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/genomic_features_1Mb/Hi-C

# get score
ln -s ../../genomic_features/Hi-C/compart_scores.100k.txt .

# make a bed with score
awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,".",$5,"."}' compart_scores.100k.txt | \
  sort-bed - > DSG.bed
awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,".",$6,"."}' compart_scores.100k.txt | \
  sort-bed - > FA.bed

# ave 1M
bedmap --echo --mean --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed DSG.bed > DSG.ave1M
bedmap --echo --mean --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed FA.bed > FA.ave1M

# paste into one
paste -d '\t' DSG.ave1M FA.ave1M | \
  cut -f 1-4,8 > compart_scores.1M.txt




