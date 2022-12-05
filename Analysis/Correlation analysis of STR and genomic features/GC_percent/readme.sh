#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/genomic_features_1Mb/GC_percent

# GC percent
bedtools nuc -fi ~/RefData/Homo_sapiens/GRCh38_no_alt/genome.fa -bed ../hg38.autosomes.bin_1M.rmCen.bed > GC_percent.txt

# to 4 cols
sed '1d' GC_percent.txt | cut -f 1-3,5 > GC_percent.count


