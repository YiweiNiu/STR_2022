#!/usr/bin/env Rscript

# run
# cd /niuyw-usb-disk/Projects/STR/str_rproj
# nohup Rscript analysis/PCA.NyuWa.str_length.R > analysis/PCA.NyuWa.str_length.log &

DOCNAME = "PCA.NyuWa.str_length"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

# Tidyverse
library(tidyverse)
# plot
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)
library(ggcorrplot)
library(ggrepel)
library(ggrastr)
library(EnvStats)
theme_set(theme_cowplot(font_size = 12, # axis.title
                        rel_small = 10/12, # axis.text
                        rel_tiny = 9/12,
                        rel_large = 12/12))
# heatmap
library(pheatmap)
# color
library(ggsci)
library(RColorBrewer)
# data.table
library(data.table)

# load dat
mat_dir = '/niuyw-usb-disk/Projects/STR/PCA/NyuWa'
chroms = paste0("chr", 1:22)
mat_lst = lapply(chroms, function(x){
  fname = file.path(mat_dir, paste0(x, ".length.csv"))
  read_csv(fname)
})
mat = do.call(rbind, mat_lst)
mat[1:3, 1:3]
rm(mat_lst)
# neat
dat = mat %>%
  column_to_rownames("Allele") %>%
  as.matrix()
dat[1:3, 1:3]
rm(mat)
# replace NA with rowMean
k <- which(is.na(dat), arr.ind=TRUE)
dat[k] <- rowMeans(dat, na.rm=TRUE)[k[,1]]
# pca
pc = prcomp(t(dat))
pc_res = data.frame(sample_old = rownames(pc$x),
                    pc1 = pc$x[,1], pc2 = pc$x[,2],
                    pc3 = pc$x[,3], pc4 = pc$x[,4],
                    pc5 = pc$x[,5], pc6 = pc$x[,6],
                    pc7 = pc$x[,7], pc8 = pc$x[,8],
                    pc9 = pc$x[,9], pc10 = pc$x[,10])
pc_var = tibble(eigs = pc$sdev^2) %>%
  mutate(pcs = 1:n(),
         Proportion = eigs/sum(eigs),
         Cumulative = cumsum(eigs)/sum(eigs))
saveRDS(pc_var, file = here::here("output", DOCNAME, "NyuWa.pc_var_df.rds"))
# save
saveRDS(pc, file = here::here("output", DOCNAME, "NyuWa.prcomp_res.rds"))
saveRDS(pc_res, file = here::here("output", DOCNAME, "NyuWa.pc_res_df.rds"))



