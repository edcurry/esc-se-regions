####################################################################
# 
# Wout Megchelenbrink
# Jan 18, 2020
#
# Enriched TF motifs vs expressed TFs in PU, DM and PM regions
#
#####################################################################

rm(list=ls())
source("include/style.R")

# The superenhancelets were tiled into regions with the size of a typical ATAC-seq peak (~300 bp). Motif analysis was done in these tiles
DT <- fread("DATA/TF_Motifs_SElets_GimmeMotifs.tsv")

DT[, tf:=tolower(tf)]


### Get the RNA=-seq expression per gene
rnaseq <- unique(fread("DATA/RNAseq.tsv")[, .(tf=tolower(gene_name), rpkm_ser=log2(rpkm_ser+1)) ])

# Merge motif TFs and RNA-seq
DTX <- merge(DT, rnaseq, by="tf",all.x=T)
DTX[is.na(rpkm_ser), rpkm_ser:=0]

# Take the TF with the lowest p-value (TFs can reside in multiple GIMME motif clusters; we take the most significant one)
DTX <- DTX[DTX[, .I[pval <= min(pval)], by=.(tf, type)]$V1,]

# Make unique to remove possible ties
DTX <- unique(DTX, by=c("tf", "type"))
setnames(DTX, "type", "SElet")

# Split TFs in to enriched only or enriched and Expressed (RPMK >= 1) 
ExpressedAndSignif <- DTX[fdr < 0.05 & rpkm_ser >= 1, .(Nexp=length(unique(tf))), by=SElet]
Signif <- DTX[fdr < 0.05, .(Nenrich=length(unique(tf))), by=SElet]

DTY <- merge(ExpressedAndSignif, Signif, by = "SElet")
DTY[, pct:=Nexp/Nenrich*100]
DTY[, pct.not:=100-pct]

# Create long format table with Nnot = enriched but not expressed and Nexp = enriched + expressed
DTY[, Nnot:=Nenrich-Nexp]
DTY.N  <- melt.data.table(DTY, id.vars = "SElet", measure.vars = c("Nexp", "Nnot"))
DTY.N[, SElet:=factor(SElet, levels = rev(c("PU", "DM", "PM")))]
DTY.N[, variable:=factor(variable, levels=c("Nnot", "Nexp"))]

# Make the barplot
ggplot(DTY.N, aes(x=SElet, y=value, fill=variable))+
geom_bar(stat = "identity", color = "#BBBBBB") +
scale_fill_manual(values = c("#D3D3D3", "#8080ff")) +
xlab("") +
ylab("") +
coord_flip() +
theme_SE() 

# Save the PDF
ggsave("IMG/FigS7A_top.pdf", width = 3.5, height=2)


#################### PIE CHARTS #################### 
DTY.pct  <- melt.data.table(DTY, id.vars = "SElet", measure.vars = c("pct", "pct.not"))
DTY.pct[, SElet:=factor(SElet, levels = c("PU", "DM", "PM"))]

# Make the pie charts
ggplot(DTY.pct, aes(x="", y=value, fill=variable))+
geom_bar(width = 1, stat = "identity", color = "#BBBBBB") +
coord_polar("y", start=0) +
scale_fill_manual(values = c("#8080ff", "#D3D3D3")) +
facet_wrap(~SElet) +
theme_void()

# Save the PDF
ggsave("IMG/FigS7A_bottom.pdf", width = 6, height = 3)
