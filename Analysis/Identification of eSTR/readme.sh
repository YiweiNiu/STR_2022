#!/bin/bash

# WORK DIR
WORKDIR=/home2/niuyw/project/STR/eSTR

# tools dir
TOOLDIR=/home2/niuyw/software
path2bcftools=$TOOLDIR/bcftools-1.14/bcftools

# get each chrom
for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
do
  input=/home2/niuyw/project/STR/02_FilteredRaw220120/${chrom}.vcf.gz
  $path2bcftools view --threads 4 -c 1 -S /home2/niuyw/project/Geuvadis/geno_pca/Geuvadis_olp_1KG2504.lst -Oz $input -o $WORKDIR/tmp.${chrom}.vcf.gz
  $TOOLDIR/htslib.1.14/bin/tabix -p vcf $WORKDIR/tmp.${chrom}.vcf.gz
done

# merge vcf
$path2bcftools concat --threads 4 -a -Oz -o $WORKDIR/Geuvadis_445.vcf.gz $WORKDIR/tmp.chr*.vcf.gz
# index
$TOOLDIR/htslib.1.14/bin/tabix -p vcf $WORKDIR/Geuvadis_445.vcf.gz

# clean
rm $WORKDIR/${pop}.chr*.vcf.gz

# statSTR
/home2/shiyr/miniconda3/bin/statSTR --vcf $WORKDIR/Geuvadis_445.vcf.gz --out $WORKDIR/Geuvadis_445 --vcftype gangstr --thresh --afreq --acount --hwep --het --entropy --mean --mode --var --numcalled

# get STRs within 500kb of genes expressed in LCL data set
awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=0}NR>FNR{if($4 in A){print $0}}' Geuvadis_445.expressed_genes.txt ~/RefData/Homo_sapiens/GENCODE_v34/gencode.v34.genes.bed6 | sort -k1,1 -k2,2n > Geuvadis_445.expressed_genes.bed6
sed '1d' Geuvadis_445.tab | cut -f 1-3 | sort -k1,1 -k2,2n | bedtools window -a Geuvadis_445.expressed_genes.bed6 -b stdin -w 500000 > Geuvadis_445.genes_window500k_str.txt

# filter STRs and normalize STR genotypes
python filter_str_for_eQTL.py --het 0.1 --numcalled 50 --tab Geuvadis_445.tab --pair Geuvadis_445.genes_window500k_str.txt --vcf Geuvadis_445.vcf.gz -o Geuvadis_445

# T01
Rscript eSTR.identification.R
Rscript eSTR.identification.permuted.R
#eSTR.identification.Rmd
#eSTR.examples.Rmd
#eSTR.cmp_GTEx.Rmd
#eSTR.cmp_Gymrek_2016.Rmd

# make final table (publication version)
awk 'BEGIN{OFS="\t"; print "gene", "site", "str.motif", "beta", "beta.se", "linreg.pval", "qval.gene", "signif.estr"}NR==FNR{FS="\t";site=$1":"$2;A[site]=$5}NR>FNR{FS=",";if($8+0<0.1){signif="True"}else{signif="False"};print $1,$2,A[$2],$3,$4,$6,$8,signif}' ~/RefData/STR/GangSTR-references/hg38_ver13.bed df.eSTR.csv | grep -v 'gene,site,'> eSTR.publish.txt


# GWAS Catalog to gmt
cd /home2/niuyw/RefData/GWAS_Catalog
python gwasCatalog_genes_to_gmt.py



