#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/diff_var

cd $WORKDIR

# sd diff
# AFR vs SAS
for i in 22
for i in {1..21}
do
  qsub get_diff_sd.pbs -N chr${i}.AFRvsSAS -v chrom=chr${i},pop1=AFR,pop2=SAS
done

# AMR vs SAS
for i in {1..22}
do
  qsub get_diff_sd.pbs -N chr${i}.AMRvsSAS -v chrom=chr${i},pop1=AMR,pop2=SAS
done

# EUR vs SAS
for i in {1..22}
do
  qsub get_diff_sd.pbs -N chr${i}.EURvsSAS -v chrom=chr${i},pop1=EUR,pop2=SAS
done

# EAS vs SAS
for i in {1..22}
do
  qsub get_diff_sd.pbs -N chr${i}.EASvsSAS -v chrom=chr${i},pop1=EAS,pop2=SAS
done

# AFR vs EAS
for i in {1..22}
do
  qsub get_diff_sd.pbs -N chr${i}.AFRvsEAS -v chrom=chr${i},pop1=AFR,pop2=EAS
done

# AMR vs EAS
for i in {1..22}
do
  qsub get_diff_sd.pbs -N chr${i}.AMRvsEAS -v chrom=chr${i},pop1=AMR,pop2=EAS
done

# EUR vs EAS
for i in {1..22}
do
  qsub get_diff_sd.pbs -N chr${i}.EURvsEAS -v chrom=chr${i},pop1=EUR,pop2=EAS
done

# AFR vs EUR
for i in {1..22}
do
  qsub get_diff_sd.pbs -N chr${i}.AFRvsEUR -v chrom=chr${i},pop1=AFR,pop2=EUR
done

# AMR vs EUR
for i in {1..22}
do
  qsub get_diff_sd.pbs -N chr${i}.AMRvsEUR -v chrom=chr${i},pop1=AMR,pop2=EUR
done

# AFR vs AMR
for i in {1..22}
do
  qsub get_diff_sd.pbs -N chr${i}.AFRvsAMR -v chrom=chr${i},pop1=AFR,pop2=AMR
done




