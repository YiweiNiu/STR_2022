#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/PCA

# get sam lst
cd $WORKDIR
mkdir NyuWa EAS world
# random 100 NyuWa x3
awk '$6!="NA" && $1=="NyuWa"' ../samples/220116_sample_info.GangSTR.txt | cut -f2 | shuf -n 100 - > samNyuWa100.lst1
awk '$6!="NA" && $1=="NyuWa"' ../samples/220116_sample_info.GangSTR.txt | cut -f2 | shuf -n 100 - > samNyuWa100.lst2
awk '$6!="NA" && $1=="NyuWa"' ../samples/220116_sample_info.GangSTR.txt | cut -f2 | shuf -n 100 - > samNyuWa100.lst3
# world, 100 samples of NyuWa + 2504 1KGP
cat ../samples/sam1KG2504.lst samNyuWa100.lst1 > samWorld.lst1
cat ../samples/sam1KG2504.lst samNyuWa100.lst2 > samWorld.lst2
cat ../samples/sam1KG2504.lst samNyuWa100.lst3 > samWorld.lst3
# EAS, 100 samples of NyuWa + 1KGP-EAS
cat <(awk '$5=="EAS" && $1=="1KGP"' ../samples/220116_sample_info.GangSTR.txt | cut -f 2) samNyuWa100.lst1 > samEAS.lst1
cat <(awk '$5=="EAS" && $1=="1KGP"' ../samples/220116_sample_info.GangSTR.txt | cut -f 2) samNyuWa100.lst2 > samEAS.lst2
cat <(awk '$5=="EAS" && $1=="1KGP"' ../samples/220116_sample_info.GangSTR.txt | cut -f 2) samNyuWa100.lst3 > samEAS.lst3
# NyuWa, 2866 samples with Province info
awk '$8!="NA" && $1=="NyuWa"' ../samples/220116_sample_info.GangSTR.txt | cut -f 2 > samNyuWa2856.lst
# NyuWa, CHN and CHS
awk '$6!="NA" && $1=="NyuWa"' ../samples/220116_sample_info.GangSTR.txt | cut -f 2 > samNyuWa.CHN_CHS.lst

# Use bcftools to subset vcfs
for i in {1..22}
do
  qsub subset_vcf.pbs -N chr${i}.subset_vcf -v chrom=chr${i}
done

# vcf to mat (length mat)
for i in {1..22}
do
  qsub vcf2mat_len.pbs -N chr${i}.vcf2mat_len -v chrom=chr${i}
done


