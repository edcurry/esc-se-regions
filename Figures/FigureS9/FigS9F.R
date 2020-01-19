######################################################################
#
# Wout Megchelenbrink
# January 18, 2020
######################################################################

rm(list=ls())
source("include/style.R")

# Parameters
pc        <- 1    # pseudo count to add 
min.size  <- 250  # only consider typical enhancers of length 250 bp or larger (for others, the small size may artificially inflate the RPKM values)

# Read and merge the data
TEs  <- fread("DATA/Typical_Enhancers.tsv")[type %in% c("PU", "DM") & End-Start >= min.size ,]
TEs[, H3K27ac_FC:=log2(H3K27ac_KO_RPKM+1)-log2(H3K27ac_FF_RPKM+1)]
TEs[, Med1_FC:=log2(Med1_KO_RPKM+1) - log2(Med1_FF_RPKM+1)]
TEs[, Esrrb:=log2(Esrrb_Exo_RPKM+1)]

# Convert to long format for ggplot
DTX <- melt.data.table(TEs, id.vars = c("GeneID", "type", "Esrrb"),
                       measure.vars = c("H3K27ac_FC", "Med1_FC"), 
                      variable.name = "chip", value.name = "log2_FC")
setnames(DTX, c("TE_id", "type", "Esrrb_log2_NC", "chip", "log2_FC"))

# Stats
unique(DTX[, cor.test(Esrrb_log2_NC, log2_FC), by=.(chip, type)][, .(chip, type, estimate, p.value)])

# Revel
DTX[, chip:=factor(chip, levels = c("Med1_FC", "H3K27ac_FC"))]
DTX[, type:=factor(type, levels = c("PU", "DM"))]

# Make the scatterplot
ggplot(DTX, aes(x=Esrrb_log2_NC, y=log2_FC, color = type, fill = type)) + #, label = r)) +
geom_point(cex = .2) +
geom_smooth(method = "lm") +
scale_color_manual(values = c("#00ff00", "#ff00ff"), guide=F) +  
scale_fill_manual(values = c("#B2F7B2", "#F28FF2"), guide=F) +  
xlab("Esrrb occupancy in WT ESC (log2 FPKM)") +
ylab("Occupancy Esrrb-/- vs WT (log2 FC)") +
scale_y_continuous(limits = c(-3,3), breaks = seq(-2,2,by=2)) +
facet_grid(~chip ~ type, scales = "free_y") +
theme_SE() 

# Save PDF
ggsave("IMG/FigS9F.pdf", width = 5, height = 5)
