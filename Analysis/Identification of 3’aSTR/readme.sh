#!/bin/bash

# T01
Rscript aSTR.identification.R
Rscript aSTR.identification.permuted.R
# aSTR.identification.Rmd
# aSTR.examples.Rmd

# BDC
cd /home2/niuyw/project/STR/aSTR

# make final table (publication version)
awk 'BEGIN{OFS="\t"; print "transcript", "gene", "site", "str.motif", "beta", "beta.se", "linreg.pval", "qval.transcript", "signif.astr"}NR==FNR{FS="\t";site=$1":"$2;A[site]=$5}NR>FNR{FS=",";if($9+0<0.1){signif="True"}else{signif="False"};print $1,$2,$3,A[$3],$4,$5,$7,$9,signif}' ~/RefData/STR/GangSTR-references/hg38_ver13.bed df.aSTR.csv | grep -v 'transcript,gene_id,'> aSTR.publish.txt

# annotate aSTRs
python anno_aSTR.py



