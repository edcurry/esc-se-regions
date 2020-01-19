######################################################################
#
# Wout Megchelenbrink
# Jan 17, 2020
# 
######################################################################

library(car)

rm(list=ls())

SElet.size.cutoff <- 500

source("include/style.R")

# Read CHIC interaction file
SElets <- fread("DATA/SElets.tsv")[type != "CG"]
SElets[ ,type:=ifelse(type == "INT", "M", "U")]
se_ids <- SElets[, .(N=length(unique(type))), by=se_id][N >= 2][, se_id]
DT <- fread("DATA/CHIC_promoter_SE_interactions_Joshi_Sahlen_629_with_RPKM_ser_above_1_and_gene_status_KNOWN.tab")[gene_type == "protein_coding" & gene_status == "KNOWN"]#[se_id %in% se_ids]
joshi.signif <- unique(DT[in.joshi == 1, .(gene_name, se_id)])
sahlen.signif <- unique(DT[in.sahlen == 1, .(gene_name, se_id)])

DT.joshi <- merge(fread("DATA/CHIC_PMT_SElet_KNOWN_genes_1RPKM_or_more.tsv"), joshi.signif, by=c("gene_name", "se_id"))[length >= SElet.size.cutoff]
DT.sahlen <- merge(fread("DATA/CHIC_PMT_SElet_KNOWN_genes_1RPKM_or_more.tsv"), sahlen.signif, by=c("gene_name", "se_id"))[length >= SElet.size.cutoff]

# Convert to log2
DT.joshi[, fpkm_2i:=log2((fpkm_2i_1+fpkm_2i_2)/2 +1)]
DT.joshi[, fpkm_ser:=log2((fpkm_ser_1+fpkm_ser_2)/2 +1)]
DT.sahlen[, fpkm_sahlen:=log2(fpkm_sahlen + 1)]

# ANOVA with blocking design on the bait count per SElet
Anova(lm(fpkm_2i ~ subtype + n_baits_joshi + log2(distance+1), data=DT.joshi), type = 3)
Anova(lm(fpkm_ser ~ subtype + n_baits_joshi + log2(distance+1) , data=DT.joshi), type =3)
Anova(lm(fpkm_sahlen ~ subtype + n_baits_sahlen + log10(distance+1) , data=DT.sahlen), type = 3)

DTX <- rbind(DT.joshi[, .(fpkm=fpkm_2i, n_baits=n_baits_joshi, subtype, class="2iL-ESC (Joshi et al.)")],
             DT.joshi[, .(fpkm=fpkm_ser, n_baits=n_baits_joshi, subtype, class="SL-ESC (Joshi et al.)")],
             DT.sahlen[, .(fpkm=fpkm_sahlen, n_baits=n_baits_sahlen, subtype, class="SL-ESC (Sahlen et al.)")])

setnames(DTX, "subtype", "SE")
DTX[, SE:=factor(SE, levels = c("PU", "DM", "INT"))]

# Create figure
ggplot(DTX, aes(x=class, y=fpkm, fill = SE)) +
geom_boxplot(width=.9, position = position_dodge(width = .9), notch = T, show.legend = T, outlier.size = .25) +
scale_fill_manual(values =c("#00ff00", "#ff00ff", "#d3d3d3")) +
xlab("Capture baits per subregion (#)") +
ylab("C-HiC intensity (log2 FPKM)") +
theme_SE() 

# Save figure
ggsave("IMG/FigS3I.pdf", width = 4, height = 3.5)
