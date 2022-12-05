
setwd('/home2/niuyw/project/STR/genomic_features/Hi-C')

library(tidyverse)
# heatmap
library(pheatmap)
library(GENOVA)

options(scipen=999)

# samples
samples = c("GSM5057481", "GSM5057489")
# load in batch
hic_100k = list(GSM5057481 = load_contacts(signal_path = 'GSM5057481_U54-ESC-DSG-HindIII-20161206-R1-T1_hg38.hg38.mapq_30.1000.mcool',
                                          sample_name = "GSM5057481",
                                          resolution = 100000),
               GSM5057489 = load_contacts(signal_path = 'GSM5057489_U54-ESC-FA-HindIII-20160311-R1-T1_hg38.hg38.mapq_30.1000.mcool',
                                          sample_name = "GSM5057489",
                                          resolution = 100000))

# DNase-seq
DNase_peaks = read.delim('/Parastor300s_G30S/niuyw/project/STR/genomic_features/DNase-seq/ENCFF905XDS.bed.gz',
                         header = FALSE)

# all
# A- and B-compartments
CS_out = compartment_score(hic_100k, bed = DNase_peaks)
saveRDS(CS_out, file = 'CS_out.100k.rds')
# CS_out = readRDS('CS_out.100k.rds')

scores = CS_out$compart_scores
write.table(scores, file = 'compart_scores.100k.txt', quote = F,
            row.names = FALSE, col.names = FALSE)


