###############################################################
# 
# Wout Megchelenbrink
# Jan 17, 2020
#
# Gene set enrichment analysis
###############################################################

rm(list=ls())
source("include/style.R")

stats <- fread("DATA/GSEA_result.tsv")[padj < 0.05]
stats[, pathway:=factor(pathway, levels=rev(stats$pathway))]

ggplot(stats, aes(x=pathway, y=score)) +
geom_bar(stat = "identity", fill = bs.col.dark[2]) + 
coord_flip() +
theme_SE()

ggsave("IMG/FigS3E.pdf", width = 8, height = 2)