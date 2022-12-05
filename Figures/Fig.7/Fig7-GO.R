library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr) 
library(ggplot2)
library(cowplot)
library(stringr)

theme_set(theme_cowplot( font_size = 20, rel_small = 10 / 10, rel_tiny = 8 / 10, rel_large = 10 / 10))

f <- read.table("gene.txt",header=T,sep=' ') 
x <-f[,1]
hg<-bitr(x,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")
go <- enrichGO(hg$ENTREZID,OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'ENTREZID')
write.csv(go,file="GO.csv")
dotplot(go,showCategory=10) + scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +theme(panel.grid=element_blank())

ggsave(filename="GO.pdf", width = 18, height = 10, units = "cm")
