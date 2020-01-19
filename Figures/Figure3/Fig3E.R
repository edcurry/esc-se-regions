######################################################################
#
# Wout Megchelenbrink
# Januari 17, 2020
# 
# ATAC-seq boxplots in ESCs, EpiLCs and EpiSCs per SElet type
######################################################################

rm(list=ls())

source("include/style.R")

# Read the ATAC-seq data
selets <- fread("DATA/SElets.tsv")[type %in% c("PU", "DM", "INT")]

# Convert to long format
DT <- melt.data.table(selets, id.vars = c("SElet_id", "type"), 
                      measure.vars =  c("ATAC_ser", "ATAC_EpiLC", "ATAC_EpiSC"),
                      variable.name = "mark")


# Annotate condition
DT[grep("ser", mark), condition:="ser-ESC"]
DT[grep("EpiLC", mark), condition:="EpiLC"]
DT[grep("EpiSC", mark), condition:="EpiSC"]
DT[, mark:=NULL]
DT[, .N, by=condition]

# Relevel factors
DT[, condition:=factor(condition, levels=c("2i-ESC", "ser-ESC", "EpiLC", "EpiSC"))]
DT[, type:=factor(type, levels=c("PU", "DM", "INT"))]

# Make the plot
ggplot(DT, aes(x=condition, y=value, fill = type)) + 
geom_boxplot(notch = T, show.legend = T, cex = .25, outlier.size = .25)  +
scale_fill_manual(values = c("#00ff00", "#ff00ff", "#d3d3d3")) +  
theme_SE() + 
xlab("") + 
ylab("") + 
theme(legend.title = element_blank())

# Save PDF
ggsave("IMG/Fig3E.pdf", width = 3.25, height = 2.5)
