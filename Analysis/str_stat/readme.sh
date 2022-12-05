#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/str_stat

# sample callrate
cd $WORKDIR/sample_callrate
for i in {1..22}
do
  qsub sample_callrate_NyuWa.pbs -N chr${i}_callrate_NyuWa -v chrom=chr${i}
  qsub sample_callrate_1KGP.pbs -N chr${i}_callrate_1KGP -v chrom=chr${i}
done
python sample_callrate_collect.py *.txt > sample_callrate.txt

# saturation curve
cd $WORKDIR/str_saturation
qsub str_saturation.pbs
qsub str_saturation.NyuWa.pbs
qsub str_saturation.1KGP.pbs

# sample stat
cd $WORKDIR/sample_stat
for i in {1..22}
do
  qsub sample_stat.pbs -N chr${i}.sample_stat -v chrom=chr${i}
done
python sample_stat_collect.py chr*.csv > sample_stat.csv

# add AC\AF\MAF 这些信息
cd /home2/niuyw/project/STR/samples
sed '1d' 220116_sample_info.GangSTR.txt | awk 'BEGIN{FS=OFS="\t"}{print $2,$1}' > sam_2_dataset.txt
# EAS not in NyuWa!
sed '1d' 220116_sample_info.GangSTR.txt | grep -v 'NyuWa' | awk 'BEGIN{FS=OFS="\t"}{print $2,$5}' > sam_2_SPop.txt
sed '1d' 220116_sample_info.GangSTR.txt | awk 'BEGIN{FS=OFS="\t"}{if($6!="NA"){print $2,$6}}' > sam_2_Pop.txt
cd $WORKDIR/AF
for i in {1..22}
do
  qsub getMoreAF.pbs -N chr${i}.getMoreAF -v chrom=chr${i}
done
python collect_AF_dataset.py chr*dataset.txt > str_af.dataset.csv
python collect_AF_spop.py chr*superPop.txt > str_af.superPop.csv
python collect_AF_pop.py chr*pop.txt > str_af.pop.csv
# get MAF (major allele frequency) of dataset and superpop
python get_maf.py

# pSTR info to table
cd $WORKDIR/str_info
for i in {1..22}
do
  qsub getMoreINFO.pbs -N chr${i}.getMoreINFO -v chrom=chr${i}
done
python collect_info.py chr*.txt > str_info.csv

# npSTR info to table
cd $WORKDIR/npstr_info
for i in {1..22}
do
  #qsub getMoreINFO.pbs -N chr${i}.getMoreINFO -v chrom=chr${i}
  bash getMoreINFO.pbs chr${i}
done
python collect_info.py chr*.txt > npstr_info.csv


# get a lst of STR in each pop
cd $WORKDIR/sample_stat_byPop
python site_lst_eachPop.py
python allele_lst_eachPop.py
# pop stat
python sample_stat_byPop.py > sample_stat_byPop.csv


# pop sharing
cd $WORKDIR/pop_sharing
python pop_sharing.py > pop_sharing.csv






