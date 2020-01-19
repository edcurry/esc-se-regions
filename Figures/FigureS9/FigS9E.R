####################################################################
# 
# Wout Megchelenbrink
# January 18, 2020
# Perform motif enrichment analysis in PM, DM and PU subregions
#
#####################################################################

rm(list=ls())
source("include/style.R")

fdr.cutoff  <- 0.05
selected    <- c("Pou5f1", "Sox2", "Klf2", "Klf4", "Klf5", "Esrrb", "Nr5a2", "Nr0b1")

# Get the RNAseq
RNAseq <- fread("DATA/RNAseq.tsv")[, .(tf=gene_name, 
                                       rpkm_ser=log2(rpkm_ser+1),
                                       rpkm_episc=log2(rpkm_episc+1))]

# Get the PU/DM typical enhancers
DT <- fread("DATA/TF_Motifs_TE_GimmeMotifs.tsv")

DTX <- merge(DT, RNAseq, by="tf", all.x=T)
DTX[is.na(rpkm_ser), rpkm_ser:=0]
DTX[is.na(rpkm_episc), rpkm_episc:=0]

# Select the most enriched occurence of each TF (TFs may be represented by multiple GIMME motifs)
DTX <- DTX[DTX[, .I[pval <= min(pval)], by=.(tf, type)]$V1,]
DTX <- unique(DTX, by=c("tf", "type"))
setnames(DTX, "type", "SElet")

DTX <- DTX[fdr < fdr.cutoff]
DTX[, score:=-log10(fdr)]
DTX[, tf:=str_to_title(tf)]
RNAseq[, tf:=str_to_title(tf)]

mat.motif.selected <- acast(DTX[tf %in% selected], formula = tf ~ SElet, value.var = "score", fill = 0)[selected, ]
mat.motif.selected <- mat.motif.selected[, c("PU", "DM")]

tmp <- RNAseq[tf %in% selected, .(tf, ESC=rpkm_ser, EpiSC=rpkm_episc)]
mat.rnaseq.selected <- data.matrix(tmp[, 2:3])
rownames(mat.rnaseq.selected) <- tmp$tf
mat.rnaseq.selected <- mat.rnaseq.selected[selected, ]

hm1 <- Heatmap(mat.motif.selected, name="Motif enrichment\n-log10(FDR)", cluster_columns = F, show_row_names = T,
               column_title = "TF motifs", rect_gp = chm.border.style.light,
               row_names_gp = gpar(fontsize = 11), row_names_side = "left",  cluster_rows = F,
               col = colorRamp2(c(0, 20), c("white", colors.dark[4])))

hm2 <- Heatmap(mat.rnaseq.selected, name="RNAseq\nlog2(RPKM)", cluster_columns = F, show_row_names = F, 
               column_title = "RNA-seq", rect_gp = chm.border.style.light,
               col = colorRamp2(c(0, 3.5, 7), c(colors.dark[3], "white", colors.dark[5])))

pdf("IMG/Fig9E.pdf", width = 4, height = 3)
hm1 + hm2
dev.off()
