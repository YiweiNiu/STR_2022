#Load packages

# Tidyverse
library(tidyverse)

# plot
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)
library(ggcorrplot)
library(ggbeeswarm)
library(ggrepel)
library(ggrastr)
library(EnvStats)
theme_set(theme_cowplot(
  font_size = 10, # Overall font size
  rel_small = 10 / 10, # axis tick labels
  rel_tiny = 8 / 10, # caption
  rel_large = 10 / 10 # title
))
geom.text.size <- 7

# heatmap
library(pheatmap)
library(corrplot)

# color
library(ggsci)
library(RColorBrewer)

#1   Load chromosome size
chrom_size = read.table("chromband-q.bed",sep=' ',header=T) %>%
  dplyr::select(chrom = Chr, start = Start, end = End)
#Load pSTR
df_pstr = read_csv("str_info_q.csv") %>%
  dplyr::select(chrom = CHROM, pos = POS, period = PERIOD, motif = RU)
#Load npSTR
df_mstr = read_csv("npstr_info_q.csv") %>%
  dplyr::select(chrom = CHROM, pos = POS, period = PERIOD, motif = RU)
# load chrom_exit
chrom_exit = read.table("chromexit-q.txt",sep=" ",header=T)

#2   Neat
#pSTR
df_pstr = df_pstr %>%
  left_join(chrom_size, by = "chrom") %>%
  mutate(dist= end-pos, 
         motif = str_to_upper(motif)) %>% 
  dplyr::select(chrom, pos, motif, period, dist)
# mSTR
df_mstr = df_mstr %>%
  left_join(chrom_size, by = "chrom") %>%
  mutate(dist=end-pos,         
         motif = str_to_upper(motif)) %>%
  dplyr::select(chrom, pos, motif, period, dist)

#3  Plot
p1 = df_pstr %>%
  mutate(bin = cut(dist, c(0:150*10^6), labels = as.character(1:150))) %>%
  group_by(period, bin) %>%
  summarise(n = n()) %>%
  mutate(chrom_count=chrom_exit[which(chrom_exit$bin_serial %in% bin),'len_chrom_exit']) %>%
  mutate(mean=n/chrom_count) 

p1= p1[which(p1$chrom_count >= 5), ]
write.table(p1,file="tele-q-pstr.txt",quote=F,col.name=T,row.names=F)

p1 = p1 %>%
  ggplot(aes(x = bin, y = mean, group = 1)) +
  geom_line() +
  scale_x_discrete(breaks = as.character(seq(0, 110, 10))) +
  labs(x = "Telomere distance (Mbp)", y = "Average number of pSTRs per 1 Mbp") +
  geom_vline(xintercept = c("5"), color = "red") +
  facet_grid(period~., scales = "free") +
  theme(strip.background = element_rect(fill = "#f3f2f1"))

p2 = df_mstr %>%
  mutate(bin = cut(dist, c(0:150*10^6), labels = as.character(1:150))) %>%
  group_by(period, bin) %>%
  summarise(n = n()) %>%
  mutate(chrom_count=chrom_exit[which(chrom_exit$bin_serial %in% bin),'len_chrom_exit']) %>%
  mutate(mean=n/chrom_count)

p2= p2[which(p2$chrom_count >= 5), ]
write.table(p2,file="tele-q-mstr.txt",quote=F,col.name=T,row.names=F)

p2 = p2 %>%
  ggplot(aes(x = bin, y = mean, group = 1)) +
  geom_line() +
  scale_x_discrete(breaks = as.character(seq(0, 110, 10))) +
  labs(x = "Telomere distance (Mbp)", y = "Average number of mSTRs per 1 Mbp") +
  geom_vline(xintercept = c("5"), color = "red") +
  facet_grid(period~., scales = "free") +
  theme(strip.background = element_rect(fill = "#f3f2f1"))

p.dist = p1 + p2
ggsave2(p.dist,
        filename = "tele-q.pdf",
        height = 6, width = 7
)
