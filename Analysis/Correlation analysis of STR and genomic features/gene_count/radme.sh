#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/genomic_features_1Mb/gene_count

ln -s ../../genomic_features/gene_count/*.bed6 .

# count number per 100-kb window
bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed all.genes.bed6 > all_gene.count
bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed protein_coding.genes.bed6 > pcg.count

