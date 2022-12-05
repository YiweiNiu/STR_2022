#!/bin/bash

# get a union set of variants
cd /home2/niuyw/RefData/STR/ExpansionHunterVariantCatalog
python union_variant_catalog.py
# then edit it manually

# to BED
python variant_catalog_json_to_bed.py # variant_catalog.bed
# sort
cat variant_catalog.bed | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > variant_catalog.sorted.bed

# cmp with GangSTR ref
cut -f 1-6 variant_catalog.sorted.bed | \
  bedtools window -a stdin -b ~/RefData/STR/GangSTR-references/hg38_ver13.bed -w 10 > tmp.txt

# edit manually，get loci overlapped with GangSTR
# variant_catalog.olp_withGangSTR.txt

# test ExpansionHunter
qsub HG00315.EH.pbs

# prepare samples for ExpansionHunter
cd /home2/niuyw/project/STR/ExpansionHunter/sh
awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($1!="sample_old"){A[$1]=$4}}NR>FNR{if($1!="dataset"){if($4=="F"){sex="female"}else{sex="male"}{print $2,sex,A[$2]}}}' ../../samples/20220113.bam.index.txt ../../samples/220116_sample_info.GangSTR.txt > sam6487.4EH

# run ExpansionHunter and str_analysis
# for 6487 samples
# NyuWa sex info was inferred
cat eh_aa | while read line
cat eh_ab | while read line
cat eh_ac | while read line
cat eh_ad | while read line
cat eh_ae | while read line
cat eh_af | while read line
cat eh_ag | while read line
cat eh_ah | while read line
do
  sample="$(cut -f1 <<< "$line")"
  sex="$(cut -f2 <<< "$line")"
  bam="$(cut -f3 <<< "$line")"
  qsub HG00315.EH.pbs -N ${sample}.EH -v sample=${sample},sex=$sex,bam=$bam -q Blade
done

# combine EH
cd /home2/niuyw/project/STR/ExpansionHunter/ExpansionHunter
~/software/anaconda3/bin/combine_expansion_hunter_json_to_tsv -c ~/RefData/STR/ExpansionHunterVariantCatalog/variant_catalog.json -c ~/RefData/STR/ExpansionHunterVariantCatalog/Stranger.variant_catalog.json --include-all-fields -o sam6487_EH
cp sam6487_* ../

# combine str-analysis
cd /home2/niuyw/project/STR/ExpansionHunter/str_analysis
~/software/anaconda3/bin/combine_json_to_tsv -d -f -o sam6487_str_analysis
cp sam6487_* ../

# cmp EH and GangSTR
cd /home2/niuyw/project/STR/ExpansionHunter/
python cmp_ExpansionHunter_gangstr.py > cmp_ExpansionHunter_gangstr.txt

# allele distribution
cd /home2/niuyw/project/STR/ExpansionHunter
python make_gnomAD_genotype.py

# get allele distribution
python get_allele_distribution.py


