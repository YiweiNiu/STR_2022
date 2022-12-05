#!/usr/bin/env Rscript

# cd /niuyw-usb-disk/Projects/STR/str_rproj
# nohup Rscript analysis/aSTR.identification.R > analysis/aSTR.identification.log &

DOCNAME <- "aSTR.identification"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

# Tidyverse
library(tidyverse)
library(broom)

# plot
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)
library(ggcorrplot)
library(ggrepel)
library(ggrastr)
library(EnvStats)
theme_set(theme_cowplot(
  font_size = 12, # axis.title
  rel_small = 10 / 12, # axis.text
  rel_tiny = 9 / 12,
  rel_large = 12 / 12
))

# heatmap
library(pheatmap)
library(corrplot)

# color
library(ggsci)
library(RColorBrewer)

# load pudi
rna_pudi <- readRDS(here::here("output", "aSTR.preprocess_pudi", "Geuvadis.pudi_adjusted.rds"))
rna_pudi[1:3, 1:3]

# load dapars2
df.DaPars2 = read_delim("/niuyw-usb-disk/Projects/Geuvadis/DaPars2/DaPars2.res.txt", delim = "\t")
df.DaPars2[1:3, 1:5]

# get gene-transcript pairs
df.gene_transcript_pair = df.DaPars2 %>%
  mutate(gene_id = str_split(Gene, "[|]", simplify = TRUE)[,2],
         transcript = str_split(Gene, "[|]", simplify = TRUE)[,1],) %>%
  filter(transcript %in% colnames(rna_pudi)) %>%
  dplyr::select(gene_id, transcript)
head(df.gene_transcript_pair)

# Load STR geno
str_geno <- read_csv("/niuyw-usb-disk/Projects/STR/eSTR/Geuvadis_445.str_geno.csv") %>%
  column_to_rownames("site")
str_geno[1:3, 1:3]

# STR-gene pairs
df.str_gene_pair <- read_csv("/niuyw-usb-disk/Projects/STR/eSTR/Geuvadis_445.str_gene_pair.csv")
head(df.str_gene_pair)

# STR-transcript pairs
df.str_transcript_pair = df.str_gene_pair %>%
  left_join(df.gene_transcript_pair, by = "gene_id")

# transcript-str lst
transcript_str_lst <- split(df.str_transcript_pair$site, df.str_transcript_pair$transcript)

# test
s <- "chr7:127498642"
t <- "ENST00000000233.10"
s_geno <- na.omit(unlist(str_geno[s, ]))
t_pudi <- rna_pudi[, t]
names(t_pudi) <- rownames(rna_pudi)
t_pudi <- t_pudi[names(s_geno)]
dat <- data.frame(sample = names(t_pudi), exp = scale(t_pudi), geno = scale(s_geno))
x <- glm(formula = exp ~ geno, data = dat)

# transcripts
transcripts <- names(transcript_str_lst)

# fit
fit_res <- lapply(1:length(transcript_str_lst), function(i) {
  t <- transcripts[i]
  fit_res_lst <- lapply(transcript_str_lst[[i]], function(s) {
    # prep
    s_geno <- na.omit(unlist(str_geno[s, ]))
    t_pudi <- rna_pudi[, t]
    names(t_pudi) <- rownames(rna_pudi)
    t_pudi <- t_pudi[names(s_geno)]
    dat <- data.frame(sample = names(t_pudi), exp = scale(t_pudi), geno = scale(s_geno))
    # glm
    fit <- glm(formula = exp ~ geno, data = dat)
    # tidy
    fit_df <- broom::tidy(fit) %>%
      filter(term == "geno")
    fit_df[1, 1] <- s
    # neat
    fit_df %>%
      mutate(transcript = t) %>%
      dplyr::select(transcript, site = term, estimate, std.error, statistic, p.value)
  })
  do.call(rbind, fit_res_lst) %>%
    mutate(q.value = p.adjust(p.value, method = "bonferroni")) %>%
    arrange(-abs(estimate / std.error), p.value)
})
names(fit_res) <- transcripts

# Save
saveRDS(fit_res, file = here::here("output", DOCNAME, "aSTR.pudi_str.glm_fit.rds"))


