#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/genomic_features_1Mb/STR_count

# count pSTR
bed_dir=/home2/niuyw/project/STR/chrom_density/STR-bed
for i in {2..6}
do
cat $bed_dir/pSTR.period_${i}.bed | sort-bed - | bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed - > pSTR.period_${i}.bin_count
done
cat $bed_dir/pSTR.bed | sort-bed - | bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed - > pSTR.bin_count

# paste into one file
paste -d '\t' pSTR.* | \
  awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$8,$12,$16,$20,$24}' > pSTR.count

# count npSTR
bed_dir=/home2/niuyw/project/STR/chrom_density/STR-bed
for i in {2..6}
do
cat $bed_dir/npSTR.period_${i}.bed | sort-bed - | bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed - > npSTR.period_${i}.bin_count
done
cat $bed_dir/npSTR.bed | sort-bed - | bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed - > npSTR.bin_count

# paste into one file
paste -d '\t' npSTR.* | \
  awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$8,$12,$16,$20,$24}' > npSTR.count




