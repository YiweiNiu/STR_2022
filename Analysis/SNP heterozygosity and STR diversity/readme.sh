#!/bin/bash

cd /home2/niuyw/project/STR/STR_diversity_SNP_het

# SNP het
# 1KGP
for i in {1..22}
do
  qsub snp_het_1KGP.pbs -N chr${i}.snp_het_1KGP -v chrom=chr${i}
done
# NyuWa
for i in {1..22} Y
do
  qsub snp_het_NyuWa.pbs -N chr${i}.snp_het_NyuWa -v chrom=chr${i}
done
# neat
python ratio.py
cat snp_het_1KGP.tsv snp_het_NyuWa.tsv > cat_snp_het.tsv
python process_in_population.py # output: population_hete.tsv

# STR diversity
for i in {1..22}
do
  qsub str_gtcheck.pbs -N chr${i}.str_gtcheck -v chrom=chr${i}
done
# neat
python norm.py -i gtcheck/chr*.discordance.txt

# gzip data
tar czvf gtcheck.tar.gz gtcheck


