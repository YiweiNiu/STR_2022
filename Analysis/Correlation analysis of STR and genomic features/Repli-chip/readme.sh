#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/genomic_features_1Mb/Repli-chip

# get score
ln -s ../../genomic_features/Repli-chip/*.hg38.bed .

# ave over 1M window
for i in ENCFF000KUF ENCFF000KUG ENCFF000KUH ENCFF000KUK ENCFF000KUL
do
  # ave
  sort-bed $i.probe.hg38.bed | \
    bedmap --echo --mean --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed - > $i.repli_ave1M
done

# ave
# H1
paste -d '\t' ENCFF000KUF.repli_ave1M ENCFF000KUG.repli_ave1M ENCFF000KUH.repli_ave1M | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8+$12;print $1,$2,$3,sum/3}' > H1.repli.1M.bdg

# H7
paste -d '\t' ENCFF000KUK.repli_ave1M ENCFF000KUL.repli_ave1M | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8;print $1,$2,$3,sum/2}' > H7.repli.1M.bdg


