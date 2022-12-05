library(UpSetR)
library(ggplot2)
df <- read.table("input.txt", header = TRUE, sep="\t")
upset(df, nsets = 6, nintersects = 100, text.scale =2, mb.ratio = c(0.5, 0.5))
