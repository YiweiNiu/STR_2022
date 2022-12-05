#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/genomic_features_1Mb/DNase-seq

# count number per 100-kb window
zcat ENCFF338KTY.bed.gz | \
  sort-bed - | \
  bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed - > ENCFF338KTY.bin_count

zcat ENCFF574LKL.bed.gz | \
  sort-bed - | \
  bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed - > ENCFF574LKL.bin_count

zcat ENCFF905XDS.bed.gz | \
  sort-bed - | \
  bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed - > ENCFF905XDS.bin_count

# paste into one file
paste -d '\t' *.bin_count | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8+$12;print $1,$2,$3,sum/3}' > DHS.count


