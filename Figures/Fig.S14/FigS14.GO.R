library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr) 
library(ggplot2)
library(cowplot)
library(stringr)

theme_set(theme_cowplot( font_size = 10, rel_small = 10 / 10, rel_tiny = 8 / 10, rel_large = 10 / 10))

f <- read.table("tri-UTR5.txt",header=T,sep=' ') 
x <-f[,1]
hg<-bitr(x,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")
go <- enrichGO(hg$ENTREZID,OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'ENTREZID')
write.csv(go,file="tri-UTR5.csv")
dotplot(go,showCategory=15) + scale_y_discrete(labels = function(x) str_wrap(x, width = 100))

ggsave(filename="tri-UTR5.pdf", width = 21, height = 15, units = "cm")
# width = 21, 25
