#!/bin/bash

# WORK DIR
WORKDIR=/home2/niuyw/project/STR/vep
cd $WORKDIR

# run vep
for i in 22
for i in {1..21}
do
  qsub vep.pbs -N chr${i}.vep -v chrom=chr${i}
done

# run vep loci
for i in {1..22}
do
  qsub vep_loci.pbs -N chr${i}.vep_loci -v chrom=chr${i}
done

