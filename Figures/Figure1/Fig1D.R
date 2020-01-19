######################################################################
#
# Wout Megchelenbrink
# Jan 19, 2020
# 
# CHiC interaction intensity at unmethylated (U) and methylated (M) regions
######################################################################

rm(list=ls())

SElet.size.cutoff <- 500

source("include/style.R")

# Read CHIC interaction file
SElets <- fread("DATA/SElets.tsv")[type != "CG"]
SElets[ ,type:=ifelse(type == "INT", "M", "U")]
se_ids <- SElets[, .(N=length(unique(type))), by=se_id][N >= 2][, se_id]
DT <- fread("DATA/CHIC_promoter_SE_interactions_Joshi_Sahlen_629_with_RPKM_ser_above_1_and_gene_status_KNOWN.tab")[gene_type == "protein_coding" & gene_status == "KNOWN"]
joshi.signif <- unique(DT[in.joshi == 1, .(gene_name, se_id)])
sahlen.signif <- unique(DT[in.sahlen == 1, .(gene_name, se_id)])

DT.joshi <- merge(fread("DATA/CHIC_PMT_SElet_KNOWN_genes_1RPKM_or_more.tsv"), joshi.signif, by=c("gene_name", "se_id"))[length >= SElet.size.cutoff]
DT.sahlen <- merge(fread("DATA/CHIC_PMT_SElet_KNOWN_genes_1RPKM_or_more.tsv"), sahlen.signif, by=c("gene_name", "se_id"))[length >= SElet.size.cutoff]

DT.joshi[, fpkm_2i:=log2((fpkm_2i_1+fpkm_2i_2)/2 +1)]
DT.joshi[, fpkm_ser:=log2((fpkm_ser_1+fpkm_ser_2)/2 +1)]
DT.sahlen[, fpkm_sahlen:=log2(fpkm_sahlen + 1)]

# ANOVA with blocking design on the bait count per SElet
summary(lm(fpkm_2i ~ type + n_baits_joshi + log2(distance+1), data=DT.joshi))
summary(lm(fpkm_ser ~ type + n_baits_joshi + log2(distance+1) , data=DT.joshi))
summary(lm(fpkm_sahlen ~ type + n_baits_sahlen + log10(distance+1) , data=DT.sahlen))


DTX <- rbind(DT.joshi[, .(fpkm=fpkm_ser, n_baits=n_baits_joshi, type, class="SL-ESC (Joshi et al.)")],
             DT.sahlen[, .(fpkm=fpkm_sahlen, n_baits=n_baits_sahlen, type, class="SL-ESC (Sahlen et al.)")])

setnames(DTX, "type", "SE")
DTX[, SE:=factor(SE, levels = c("U", "M"))]

# Create figure
ggplot(DTX, aes(x=class, y=fpkm, fill = SE)) +
geom_boxplot(width=.9, position = position_dodge(width = .9), notch = T, show.legend = T) +
scale_fill_manual(values = c("#8080ff", "#d3d3d3")) +
xlab("Capture baits per subregion (#)") +
ylab("C-HiC intensity (log2 FPKM)") +
theme_SE() 

# Save figure
ggsave("IMG/Fig1D.pdf", width = 4, height = 3.5)
