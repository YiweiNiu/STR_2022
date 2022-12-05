#!/usr/bin/env Rscript

# cd /niuyw-usb-disk/Projects/STR/str_rproj
# nohup Rscript analysis/eSTR.identification.R > analysis/eSTR.identification.log &

DOCNAME <- "eSTR.identification"
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

# load exp
rna_exp <- readRDS(here::here("output", "eSTR.preprocess_exp", "Geuvadis.exp_adjusted.rds"))
rna_exp[1:3, 1:3]

# Load STR geno
str_geno <- read_csv("/niuyw-usb-disk/Projects/STR/eSTR/Geuvadis_445.str_geno.csv") %>%
  column_to_rownames("site")
str_geno[1:3, 1:3]

# STR-gene pairs
df.str_gene_pair <- read_csv("/niuyw-usb-disk/Projects/STR/eSTR/Geuvadis_445.str_gene_pair.csv")
head(df.str_gene_pair)

# gene-str lst
gene_str_lst <- split(df.str_gene_pair$site, df.str_gene_pair$gene_id)

# test
s <- "chr20:50852367"
g <- "ENSG00000000419.12"
s_geno <- na.omit(unlist(str_geno[s, ]))
g_exp <- rna_exp[, g]
names(g_exp) <- rownames(rna_exp)
g_exp <- g_exp[names(s_geno)]
dat <- data.frame(sample = names(g_exp), exp = scale(g_exp), geno = scale(s_geno))
x <- glm(formula = exp ~ geno, data = dat)

# genes
genes <- names(gene_str_lst)

# fit
fit_res <- lapply(1:length(gene_str_lst), function(i) {
  g <- genes[i]
  fit_res_lst <- lapply(gene_str_lst[[i]], function(s) {
    # prep
    s_geno <- na.omit(unlist(str_geno[s, ]))
    g_exp <- rna_exp[, g]
    names(g_exp) <- rownames(rna_exp)
    g_exp <- g_exp[names(s_geno)]
    dat <- data.frame(sample = names(g_exp), exp = scale(g_exp), geno = scale(s_geno))
    # glm
    fit <- glm(formula = exp ~ geno, data = dat)
    # tidy
    fit_df <- broom::tidy(fit) %>%
      filter(term == "geno")
    fit_df[1, 1] <- s
    # neat
    fit_df %>%
      mutate(gene = g) %>%
      dplyr::select(gene, site = term, estimate, std.error, statistic, p.value)
  })
  do.call(rbind, fit_res_lst) %>%
    mutate(q.value = p.adjust(p.value, method = "bonferroni")) %>%
    arrange(-abs(estimate / std.error), p.value)
})
names(fit_res) <- genes

# Save
saveRDS(fit_res, file = here::here("output", DOCNAME, "eSTR.exp_str.glm_fit.rds"))

