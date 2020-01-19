######################################################################
#
# Wout Megchelenbrink
# Jan 18, 2020
######################################################################

rm(list=ls())
source("include/style.R")

# pseudo count
pc <- 1

# Read ChIP-seq data
selets <- fread("DATA/SElets.tsv")[type %in% c("PU","DM") & nCG > 0, ]

# Get logFC
selets[, H3K27ac_FC:=log2(H3K27ac_KO+pc)-log2(H3K27ac_WT+pc)]
selets[, Med1_FC:=log2(Med1_KO+pc)-log2(Med1_WT+pc)]
selets[, Esrrb:=log2(Esrrb_ser+pc)]

# Melt data table to long format
DTX <- melt.data.table(selets, id.vars = c("se_id","SElet_id", "type", "Esrrb"),
                       measure.vars = c("H3K27ac_FC", "Med1_FC"), 
                      variable.name = "chip", value.name = "FC")
setnames(DTX, c("SE_id", "SElet_id", "type", "Esrrb_log2_NC", "chip", "log2_FC"))


# Statistics
unique(DTX[, cor.test(Esrrb_log2_NC, log2_FC), by=.(chip, type)][, .(chip, type, estimate, p.value)])

# Relevel
DTX[, chip:=factor(chip, levels = c("Med1_FC", "H3K27ac_FC"))]
DTX[, type:=factor(type, levels = c("PU", "DM"))]

# Make the plot
ggplot(DTX, aes(x=Esrrb_log2_NC, y=log2_FC, color = type, fill = type)) + 
geom_point(cex = .2) +
geom_smooth(method = "lm") +
scale_color_manual(values = c("#00ff00", "#ff00ff"), guide=F) +  
scale_fill_manual(values = c("#B2F7B2", "#F28FF2"), guide=F) +  
xlab("Esrrb occupancy in WT ESC (log2 FPKM)") +
ylab("Occupancy Esrrb-/- vs WT (log2 FC)") +
facet_grid(~chip ~ type, scales = "free_y") +
theme_SE() 

# Save PDF
ggsave("IMG/Fig6E.pdf", width = 5, height = 5)
