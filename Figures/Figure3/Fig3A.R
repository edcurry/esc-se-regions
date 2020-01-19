################################################################################
# Wout Megchelenbrink
# Jan 19, 2019
#
# Predicted vs measured RNA expression change
################################################################################

library(ggrepel)
rm(list=ls())
source("include/style.R")

# Get the Capture Hi-C per transcript-SElet combination up to a distance of 1 MB
chic <- fread("DATA/CHIC_tx_SElet_Dpn2.tsv")[distance <= 1e6]
chic[, log2dist:=log2(distance+1)]

# Get the super enhancelets classified as PU or DM that have at least one CpG dinucleotide
SElets <- fread("DATA/SElets.tsv")[type %in% c("PU", "DM") & nCG > 0]
DT <- merge(chic, SElets, by=c("se_id", "SElet_id"))


# Get the geometric mean of the replicates
DT[, reads_ser:=sqrt(reads_ser_1*reads_ser_2)]


# Take strongest interaction per SElet - gene combination
DT <- DT[DT[, .I[reads_ser >= max(reads_ser)], by=.(SElet_id, gene_name)]$V1,]
DT <- unique(DT, by=c("SElet_id", "gene_name"))

# RNAseq
RNAseq <- fread("DATA/RNAseq.tsv")
DT <- merge(DT, RNAseq, by="gene_name")

DT <- unique(DT, by=c("SElet_id", "gene_name"))

# RNAseq in vitro
DT.in.vitro <- DT[, .(gene_name, gene_type,
                      x=log2(reads_ser+1) * (pct_meth_ser_marks-pct_meth_episc_veillard), 
                      y=log2(rpkm_episc+1) - log2(rpkm_ser+1), assay="In vitro (ESC to Epi)")]

# RNAseq in vivo
DT.in.vivo <- DT[, .(gene_name, gene_type,
                      x=log2(reads_ser+1) * (E35ICM_pct_meth-E65Epi_pct_meth), 
                      y=log2(rpkm_E65Epi+1) - log2(rpkm_E35ICM+1), assay="In vivo (E3.5 ICM to E6.5 Epi)")]

DT <- rbind(DT.in.vitro, DT.in.vivo)
DTX <- DT[, lapply(.SD, mean), by=.(gene_name, gene_type, assay), .SDcols = c("x","y")]

r <- DTX[!is.na(x) & !is.na(y), .(estimate=cor(x,y)), by=assay]

r
# r = 0.38 here (0.37 in the manuscript), because we simplified the input files, which may have rounded numbers.

# Naive and primed genes to highlight in the plot
genes.naive <- c("Prdm14", "Klf4", "Tbx3", "Esrrb", "Tcfcp2l1", "Dppa5a", "Zfp42", "Tdh")
genes.primed <- c("Otx2", "Lefty1", "Pou5f1", "Sox2", "Tet1", "Smarcad1", "Nanog", "Klf13", "Med13l")

DTX[, col:="#CCCCCC"]
DTX[, size:=.5]
DTX[gene_name %in% genes.naive, `:=` (col=alpha("#ff00ff", alpha = .55), size = 2)]
DTX[gene_name %in% genes.primed, `:=` (col=alpha("#00a800", alpha = .55), size = 2)]

DTX[col != "#CCCCCC", label:=gene_name]
source("include/style.R")

ggplot(DTX, aes(x=x, y=y, label=label, color=col, size=size)) +
geom_point(data = DTX[col =="#CCCCCC"], aes(x=x, y=y, label=label, color=col, size=size)) +
geom_point(data = DTX[col !="#CCCCCC"], aes(x=x, y=y, label=label, color=col, size=size)) +
scale_fill_manual(values = c(alpha(bs.col.dark[1], .95), alpha(bs.col.dark[2], .5)), guide=F) +
geom_hline(yintercept = 0, linetype = "dashed", color = "#555555") +
geom_vline(xintercept = 0, linetype = "dashed", color = "#555555") +
xlab("Predicted expression change (log2 FC)") +
ylab("Measured expression change (log2 FC)") +
scale_color_identity()  +
geom_text_repel(data = DTX[col == alpha("#00a800", alpha = .55)], size=4, show.legend = F, color = "#00a800") +
geom_text_repel(data = DTX[col == alpha("#ff00ff", alpha = .55)], size=4, show.legend = F, color = "#ff00ff") +
theme_SE() +
facet_wrap(~assay) +
annotate(x=-6.5, y=4.25, geom="text", label=sprintf("r=%2.2f", r$estimate), size=5) + 
scale_x_continuous(limits = c(-8, .5)) + 
scale_size_continuous(range = c(.5,2), guide = F)


ggsave("IMG/Fig3A.pdf", width = 6, height = 4)
