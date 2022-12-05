#!/usr/bin/env python

'''
purpose: annotate eSTRs with other annotations

usage:
cd /home2/niuyw/project/STR/eSTR
python anno_eSTR.py

output
df.eSTR.anno.tsv
'''

import sys
from itertools import islice
from collections import defaultdict
import logging


#VERSION = '220412' # first run
VERSION = '220510' # fix a bug to compute num_alleles

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main():
    """The main function
    """
    # ensg to gene symbol
    gene_anno = defaultdict(lambda: 'NA')
    with open('/home2/niuyw/RefData/Homo_sapiens/GENCODE_v34/anno.table', 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            gene_id, gene_symbol = line[0], line[2]
            gene_id = gene_id.split('.')[0]
            gene_anno[gene_id] = gene_symbol
    # str maf (major allele freq)
    str_maf = defaultdict(lambda: 'NA')
    with open('/home2/niuyw/project/STR/str_stat/AF/str_af.dataset.csv', 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            site = line[0] + ':' + line[1]
            af_lst = [float(i) for i in line[4].split(',')]
            maf = max([1-sum(af_lst), max(af_lst)])
            str_maf[site] = str(maf)
    # str anno
    str_anno = defaultdict(lambda: ['NA', 'NA', 'NA'])
    with open('/Parastor300s_G30S/shiyr/gangSTR/process-all/10-loci_info/loci_info_total.txt', 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split('\t')
            site = line[0] + ':' + line[1]
            num_alleles = str(len(line[4].split(',')))
            motif = line[13]
            csq, gene_id, extra = line[18:21]
            values_dict = {}
            for attr in extra.split(';'):
                attr, _, val = attr.strip().partition('=')
                values_dict[attr] = val
            if "SYMBOL" in values_dict:
                gene_symbol = values_dict["SYMBOL"]
            else:
                gene_symbol = 'NA'
            str_anno[site] = [num_alleles, motif, csq, gene_id, gene_symbol]
    # omim_genes
    omim_genes = defaultdict(lambda: 'NA')
    with open('/home2/niuyw/RefData/OMIM/rawdata/20180108/morbidmap_new.txt', 'rt') as fin:
        for line in fin:
            line = line.strip().split('\t')
            gene_id, phe_omim_id = line[1], line[3]
            if gene_id not in omim_genes:
                omim_genes[gene_id] = phe_omim_id
            else:
                omim_genes[gene_id] = omim_genes[gene_id] + ';' + phe_omim_id
    # disease loci
    str_disease = defaultdict(lambda: 'NA')
    with open('/home2/niuyw/RefData/STR/ExpansionHunterVariantCatalog/variant_catalog.olp_withGangSTR.txt', 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split('\t')
            site = line[6] + ":" + line[7]
            str_disease[site] = line[3] # disease
    # gwas risk
    str_gwas_snp = defaultdict(lambda: ['NA', 'NA'])
    with open('/home2/niuyw/project/STR/GWAS_risk/STR_anno_byGWASCatalog.txt', 'rt') as fin:
        for line in fin:
            line = line.strip().split()
            site, snp, trait = line[0], line[1], line[3]
            if site not in str_gwas_snp:
                str_gwas_snp[site] = [snp, trait]
            else:
                str_gwas_snp[site][0] = str_gwas_snp[site][0] + ';' + snp
                str_gwas_snp[site][1] = str_gwas_snp[site][1] + ';' + trait
    # in previous reports?
    cmp_GTEx = defaultdict(lambda: 'False')
    with open('/home2/niuyw/project/STR/eSTR/cmp_GTEx/OLP.ALL.txt', 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            site, gene_id, signif = line[0], line[1], eval(line[4])
            if signif == True:
                cmp_GTEx[(site, gene_id)] = str(signif)
    cmp_Gymrek = defaultdict(lambda: 'False')
    with open('/home2/niuyw/project/STR/eSTR/cmp_Gymrek_2016/OLP.ALL.txt', 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split()
            site, gene_id, signif = line[0], line[1], eval(line[4])
            if signif == True:
                cmp_Gymrek[(site, gene_id)] = str(signif)
    # output
    fout = open('df.eSTR.anno.tsv', 'wt')
    fout.write('\t'.join(["estr.gene_id", "estr.gene_symbol", "estr.gene.omim", "beta", "qval.gene",
                          "site", "motif", "num.alleles", "maf",
                          "csq", "gene_id", "gene_symbol",
                          "omim",
                          "disease",
                          "gwas.snp", "gwas.trait",
                          "in_Gymrek", "in_GTEx"]) + '\n')
    # load eSTRs
    estr = {}
    with open('df.eSTR.csv', 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split(',')
            gene_id, site, beta, p = line[0], line[1], line[2], float(line[7])
            gene_id = gene_id.split('.')[0]
            gene_symbol = gene_anno[gene_id]
            # keep only eSTRs
            if p >= 0.1:
                continue
            # gene in omim?
            gene_omim = omim_genes[gene_id]
            # str_maf
            maf = str_maf[site]
            # str_anno
            num_alleles, motif, csq, str_gene_id, str_gene_symbol = str_anno[site]
            # str in omim?
            str_omim = omim_genes[str_gene_id]
            # is this a disease STR?
            disease = str_disease[site]
            # is this site in gwas risk SNP LD?
            snp, trait = str_gwas_snp[site]
            # in previous reports?
            in_GTEx = cmp_GTEx[(site, gene_id)]
            in_Gymrek = cmp_Gymrek[(site, gene_id)]
            # out
            to_w = [gene_id, gene_symbol, gene_omim, beta, str(p),
                    site, motif, num_alleles, maf,
                    csq, str_gene_id, str_gene_symbol, str_omim,
                    disease, snp, trait,
                    in_Gymrek, in_GTEx]
            fout.write('\t'.join(to_w) + '\n')
    # close
    fout.close()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)

