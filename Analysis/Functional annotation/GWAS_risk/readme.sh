#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/GWAS_risk

# get GWAS SNPs with p <= 5E-8
python parse_gwas_catalog.py > gwas_snps.txt

# get window 250k STRs
awk 'BEGIN{FS=OFS="\t"}{print $1,$2-1,$2}' gwas_snps.txt | sort -k1,1 -k2,2n > gwas_snps.bed
sed '1d' ../str_stat/str_info/str_info.csv | awk 'BEGIN{FS=",";OFS="\t"}{print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools window -a gwas_snps.bed -b stdin -w 250000 > gwas_snps.window250k_pSTRs.txt

# get SNP genotypes
for i in `seq 1 22`
do
  bcftools view --threads 4 -R gwas_snps.nochr.bed -S ~/project/STR/samples/sam6487.lst -c 1 /Parastor300s_G30S/lius/Haplotype/haplotype_merge/chr${i}_merge_own_1kgp_hgdp_siteandhap.vcf.gz -Ov > SNPs/chr${i}.vcf
  bgzip -@ 4 SNPs/chr${i}.vcf
  tabix -p vcf SNPs/chr${i}.vcf.gz
  plink --vcf SNPs/chr${i}.vcf.gz --keep-allele-order --double-id --recode A --make-bed --out SNPs/chr${i}
done
#bcftools concat --threads 4 -a -Ov SNPs/*.vcf.gz > SNPs/gwas_snps.vcf
#plink --vcf SNPs/gwas_snps.vcf --keep-allele-order --double-id --recode A --make-bed --out SNPs/tmp

# STR to matrix
for i in `seq 1 22`
do
  #python str_to_mat.py -v ../02_FilteredRaw220120/chr${i}.vcf.gz -o STRs/chr${i}
  qsub str_to_mat.pbs -N chr${i}.str_to_mat -v chrom=chr${i}
done

# compute LDs of SNP-STR
#time python get_SNP_STR_LD.v2.py --pair gwas_snps.window250k_pSTRs.txt --chrom chr21 --snp SNPs/chr21.raw --str STRs/chr21.str_dosage.csv --output chr21
for i in `seq 1 20`
do
  #time python get_SNP_STR_LD.py --pair gwas_snps.window250k_pSTRs.txt --chrom chr${i} --snp SNPs/chr${i}.raw --str STRs/chr${i}.str_dosage.csv --output chr${i}
  qsub get_SNP_STR_LD.pbs -N chr${i}.LD -v chrom=chr${i}
done


# neat results
python neat_results.py
# STR_anno_byGWASCatalog.txt


