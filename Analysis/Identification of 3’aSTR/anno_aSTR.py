#!/usr/bin/env python

'''
purpose: annotate eSTRs with other annotations

usage:
cd /home2/niuyw/project/STR/aSTR
python anno_aSTR.py

output
df.aSTR.anno.tsv
'''

import sys
from itertools import islice
from collections import defaultdict
import logging


VERSION = '220511' # first run

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
    # str info
    str_anno = {}
    with open('/home2/niuyw/project/STR/str_info/str_info.txt', 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip('\n').split('\t')
            site, motif, num_alleles, maf = line[0], line[4], line[7], line[8]
            csq, gene_id, gene_symbol, gene_type, gene_in_omim, disease, gwas_snp, gwas_trait, eSTR, _, _, = line[15:]
            str_anno[site] = [motif, num_alleles, maf,
                              csq, gene_id, gene_symbol, gene_type,
                              gene_in_omim, disease, gwas_snp, gwas_trait, eSTR]
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
    # output
    fout = open('df.aSTR.anno.tsv', 'wt')
    fout.write('\t'.join(["astr.trans_id", "astr.gene_id", "astr.gene_symbol", "astr.gene.omim", "beta", "qval.trans",
                          "site", "motif", "num.alleles", "maf",
                          "csq", "gene_id", "gene_symbol", "gene_type",
                          "omim",
                          "disease",
                          "gwas.snp", "gwas.trait",
                          "eSTR.gene"]) + '\n')
    # load aSTRs
    with open('df.aSTR.csv', 'rt') as fin:
        for line in islice(fin, 1, None):
            line = line.strip().split(',')
            trans_id, gene_id, site, beta, p = line[0], line[1], line[2], line[3], float(line[8])
            gene_id = gene_id.split('.')[0]
            gene_symbol = gene_anno[gene_id]
            # keep only aSTRs
            if p >= 0.1:
                continue
            # gene in omim?
            gene_omim = omim_genes[gene_id]
            # str_anno
            motif, num_alleles, maf, csq, str_gene_id, str_gene_symbol, str_gene_type, str_omim, disease, snp, trait, eSTR = str_anno[site]
            # out
            to_w = [trans_id, gene_id, gene_symbol, gene_omim, beta, str(p),
                    site, motif, num_alleles, maf,
                    csq, str_gene_id, str_gene_symbol, str_gene_type, str_omim,
                    disease, snp, trait,
                    eSTR]
            fout.write('\t'.join(to_w) + '\n')
    # close
    fout.close()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)

