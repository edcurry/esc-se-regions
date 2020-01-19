######################################################################
#
# Wout Megchelenbrink
# Jan 17, 2020
######################################################################

rm(list=ls())
source("include/style.R")

# Read ATAC-seq coverage over input per SElet
SElets <- fread("DATA/SElets.tsv")[type != "CG"]

# Convert to long for ggplot
DT <- melt.data.table(SElets, id.vars = c("SElet_id", "type"), 
                      measure.vars =  c("ATAC_2i", "ATAC_ser", "ATAC_EpiSC"),
                      variable.name = "mark")

# Annotate the ATAC-seq condition
DT[grep("2i", mark), condition:="2i-ESC"]
DT[grep("ser", mark), condition:="ser-ESC"]
DT[grep("EpiSC", mark), condition:="EpiSC"]

# Make the plot
DT[, condition:=factor(condition, levels=c("2i-ESC", "ser-ESC", "EpiSC"))]
DT[, type:=factor(type, levels=c("PU", "DM", "INT"))]
ggplot(DT, aes(x=condition, y=value, fill = type)) + 
geom_boxplot(notch = T, show.legend = T, cex = .25, outlier.size = .25)  +
scale_fill_manual(values = c("#00ff00", "#ff00ff", "#d3d3d3")) +  
theme_SE() + 
xlab("") + 
ylab("ATAC (log2 ATAC/input)") + 
theme(legend.title = element_blank())

# Save PDF
ggsave("IMG/FigS3J.pdf", width = 3.25, height = 2.5)
