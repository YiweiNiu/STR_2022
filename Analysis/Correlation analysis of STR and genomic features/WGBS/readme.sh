#!/bin/bash

# get methy level
Rscript methylKit.R

# autosomes
awk 'BEGIN{FS=OFS="\t"}{if(/^chr[0-9]+\t/){print $1,$2-1,$3,$4}}' H1.CpG.bdg | \
  sort -k1,1 -k2,2n > H1.CpG.autosomes.bdg

awk 'BEGIN{FS=OFS="\t"}{if(/^chr[0-9]+\t/){print $1,$2-1,$3,$4}}' H1.CHG.bdg | \
  sort -k1,1 -k2,2n > H1.CHG.autosomes.bdg

awk 'BEGIN{FS=OFS="\t"}{if(/^chr[0-9]+\t/){print $1,$2-1,$3,$4}}' H1.CHH.bdg | \
  sort -k1,1 -k2,2n > H1.CHH.autosomes.bdg



