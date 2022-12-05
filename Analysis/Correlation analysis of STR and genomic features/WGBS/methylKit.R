#!/usr/bin/env Rscript

library(methylKit)
#Sys.setenv('R_MAX_VSIZE'=100000000000)
library(tidyverse)

files <- list("ENCFF601NBW.bed.gz", "ENCFF918PML.bed.gz")
cpg <- methRead(files,
            sample.id = list("ENCFF601NBW", "ENCFF918PML"),
            assembly = 'hg38',
            treatment = c(0, 0),
            header=FALSE,
            context = 'CpG',
            mincov = 5,
            resolution = 'base',
            pipeline = list(fraction = FALSE,
                            chr.col = 1,
                            start.col = 3,
                            end.col = 3, 
                            coverage.col = 5,
                            strand.col = 6,
                            freqC.col = 11
            ))
saveRDS(cpg, file = 'CpG.rds')

# 100K
tiles_cpg = tileMethylCounts(cpg, win.size=100000, step.size=100000, cov.bases = 500, mc.cores = 4)
meth_cpg = methylKit::unite(tiles_cpg)
pool_cpg = pool(meth_cpg, sample.ids=c("0"))
data.frame(chr=pool_cpg$chr, start=pool_cpg$start, end=pool_cpg$end, perc.meth=pool_cpg$numCs1/pool_cpg$coverage1) %>%
  write_tsv(file = 'H1.CpG.bdg', col_names = FALSE)


files <- list("ENCFF379ZXG.bed.gz", "ENCFF524BMX.bed.gz")
chg <- methRead(files,
            sample.id = list("ENCFF379ZXG", "ENCFF524BMX"),
            assembly = 'hg38',
            treatment = c(0, 0),
            header=FALSE,
            context = 'CHG',
            resolution = 'base',
            mincov = 5,
            pipeline = list(fraction = FALSE,
                            chr.col = 1,
                            start.col = 3,
                            end.col = 3, 
                            coverage.col = 5,
                            strand.col = 6,
                            freqC.col = 11
            ))
tiles_chg = tileMethylCounts(chg, win.size=100000, step.size=100000, cov.bases = 2500, mc.cores = 4)
meth_chg = methylKit::unite(tiles_chg)
pool_chg = pool(meth_chg, sample.ids=c("0"))
data.frame(chr=pool_chg$chr, start=pool_chg$start, end=pool_chg$end, perc.meth=pool_chg$numCs1/pool_chg$coverage1) %>%
  write_tsv(file = 'H1.CHG.bdg', col_names = FALSE)

files <- list("ENCFF417VRB.bed.gz", "ENCFF086MMC.bed.gz")
chh <- methRead(files,
            sample.id = list("ENCFF417VRB", "ENCFF086MMC"),
            assembly = 'hg38',
            treatment = c(0, 0),
            header=FALSE,
            context = 'CHH',
            mincov = 5,
            resolution = 'base',
            pipeline = list(fraction = FALSE,
                            chr.col = 1,
                            start.col = 3,
                            end.col = 3, 
                            coverage.col = 5,
                            strand.col = 6,
                            freqC.col = 11
            ))
saveRDS(chh, file = 'CHH.rds')
tiles_chh = tileMethylCounts(chh, win.size=100000, step.size=100000, cov.bases = 10000, mc.cores = 4)
meth_chh = methylKit::unite(tiles_chh)
pool_chh = pool(meth_chh, sample.ids=c("0"))
data.frame(chr=pool_chh$chr, start=pool_chh$start, end=pool_chh$end, perc.meth=pool_chh$numCs1/pool_chh$coverage1) %>%
  write_tsv(file = 'H1.CHH.bdg', col_names = FALSE)


