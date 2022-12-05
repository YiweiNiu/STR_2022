#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/genomic_features_1Mb

# bin autosomes into 1Mb windows
awk 'BEGIN{FS=OFS="\t"}/^chr[0-9]+\t/{print $1, $2}' ~/RefData/Homo_sapiens/hg38.chrom.sizes > hg38.autosomes.sizes
bedtools makewindows -g hg38.autosomes.sizes -w 1000000 > hg38.autosomes.bin_1M.bed

# remove centromere overlaps
bedtools intersect -a hg38.autosomes.bin_1M.bed -b ~/RefData/SV_filter_regions/centromeres.hg38.bed -v | sort-bed - > hg38.autosomes.bin_1M.rmCen.bed

# Find the script in each sub-dir

