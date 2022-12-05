#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/fold_change

# split samples
cd /home2/niuyw/project/STR/samples
# No NyuWa in EAS
for pop in AFR AMR EAS EUR SAS
do
  awk -F "\t" -v name="$pop" '($2 == name){print $1}' sam_2_SPop.txt > lst.${pop}
done

cd $WORKDIR

# AFR vs SAS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.AFRvsSAS -v chrom=chr${i},pop1=AFR,pop2=SAS
done

# AMR vs SAS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.AMRvsSAS -v chrom=chr${i},pop1=AMR,pop2=SAS
done

# EUR vs SAS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.EURvsSAS -v chrom=chr${i},pop1=EUR,pop2=SAS
done

# EAS vs SAS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.EASvsSAS -v chrom=chr${i},pop1=EAS,pop2=SAS
done

# AFR vs EAS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.AFRvsEAS -v chrom=chr${i},pop1=AFR,pop2=EAS
done

# AMR vs EAS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.AMRvsEAS -v chrom=chr${i},pop1=AMR,pop2=EAS
done

# EUR vs EAS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.EURvsEAS -v chrom=chr${i},pop1=EUR,pop2=EAS
done

# AFR vs EUR
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.AFRvsEUR -v chrom=chr${i},pop1=AFR,pop2=EUR
done

# AMR vs EUR
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.AMRvsEUR -v chrom=chr${i},pop1=AMR,pop2=EUR
done

# AFR vs AMR
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.AFRvsAMR -v chrom=chr${i},pop1=AFR,pop2=AMR
done

# CHNvsCHS
cd $WORKDIR
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.CHNvsCHS -v chrom=chr${i},pop1=CHN.NyuWa,pop2=CHS.NyuWa
done

# EAS
# CHB vs CHS
cd $WORKDIR
for i in {1..22}
do
  bash get_STR_FC.sh chr${i} CHB.1KGP CHS.1KGP
done
# CDX vs CHS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.CDXvsCHS.1KGP -v chrom=chr${i},pop1=CDX,pop2=CHS.1KGP
done
# KHV vs CHS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.KHVvsCHS.1KGP -v chrom=chr${i},pop1=KHV,pop2=CHS.1KGP
done
# JPT vs CHS
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.JPTvsCHS.1KGP -v chrom=chr${i},pop1=JPT,pop2=CHS.1KGP
done
# CHB vs JPT
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.CHB.1KGPvsJPT -v chrom=chr${i},pop1=CHB.1KGP,pop2=JPT
done
# CDX vs JPT
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.CDXvsJPT -v chrom=chr${i},pop1=CDX,pop2=JPT
done
# KHV vs JPT
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.KHVvsJPT -v chrom=chr${i},pop1=KHV,pop2=JPT
done
# CHB vs KHV
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.CHB.1KGPvsKHV -v chrom=chr${i},pop1=CHB.1KGP,pop2=KHV
done
# CDX vs KHV
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.CDXvsKHV -v chrom=chr${i},pop1=CDX,pop2=KHV
done
# CHB vs CDX
for i in {1..22}
do
  qsub get_STR_FC.pbs -N chr${i}.CHB.1KGPvsCDX -v chrom=chr${i},pop1=CHB.1KGP,pop2=CDX
done

# get STR length (for ploting)
for i in {1..21}
do
  qsub get_STR_len.pbs -N chr${i}.str_len -v chrom=chr${i}
done
gzip STR_len_mat/*.csv

