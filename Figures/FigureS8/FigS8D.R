####################################################################
# 
# Wout Megchelenbrink
# Jan 18, 2020
#
# Plot Pol2, GRO-seq and START-seq at PU, DM, PM regions in ser/2i ESCs
#
#####################################################################

rm(list=ls())
DT <- fread("DATA/SElets.tsv")[type %in% c("PU","DM","INT")]

DTY <- melt.data.table(DT, id.vars = c("SElet_id", "type"), measure.vars = c("GROseq_2i", "GROseq_ser", "Pol2_2i", "Pol2_ser", "STARTseq_2i", "STARTseq_ser"),
                       variable.name = "mark", value.name = "log2_rpkm")

DTY[, condition:=ifelse(str_detect(mark, "2i"), "2i", "ser")]
DTY[, condition:=factor(condition, levels=c("ser", "2i"))]
DTY[, type:=factor(type, levels=c("PU","DM", "INT"))]
DTY[, mark:=str_sub(mark, end=str_locate(mark, "_")[,1]-1)]
DTY[, mark:=factor(mark, levels=c("Pol2","GROseq", "STARTseq"))]


source("include/style.R")
ggplot(DTY, aes(x=mark, y=log2_rpkm, fill=type)) +
geom_boxplot(notch=F, cex = .25, outlier.size = .25) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
scale_fill_manual(values = c("#00ff00","#ff00ff", "#d3d3d3")) +
theme_SE() +
theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
xlab(NULL) +
facet_wrap(~condition) +
ylab("Tag density (log2 RPKM)")

ggsave("IMG/FigS8D.pdf", width = 5, height = 3)

