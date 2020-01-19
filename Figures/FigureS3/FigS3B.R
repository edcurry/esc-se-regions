##################################################################################################
# Wout Megchelenbrink
#
# Jan 18, 2020
# Scatterplot of CHIC transcript-selet interaction strength in different datasets
# 
##################################################################################################

rm(list=ls())

source("include/style.R")

# Read transcript - selet interactions
tx.selet <- fread("DATA/CHIC_tx_SElet_Dpn2.tsv")

# Convert to long format
tx.selet <- melt.data.table(tx.selet, id.vars = c("gene_name", "transcript_id",  "se_id", "SElet_id", "nb_pmt_joshi", "nb_pmt_sahlen"), 
                measure.vars = c("reads_2i_1", "reads_2i_2", "reads_ser_1", "reads_ser_2", "reads_sahlen"), 
                variable.name = "condition", value.name = "reads") 

# Merge and log2 transform
DT <- merge(tx.selet, tx.selet, by= c("gene_name", "transcript_id", "se_id", "SElet_id"), suffixes = c("_1", "_2"), allow.cartesian = T)
DT[, condition_1:=str_replace(condition_1, "reads_", "")]
DT[, condition_2:=str_replace(condition_2, "reads_", "")]
DT[, reads_1:=log2(reads_1+1)]
DT[, reads_2:=log2(reads_2+1)]

# Dicard the tx-selet interactions between Joshi and Sahlen that were not baited in both data sets
DT <- DT[!(condition_1 %in% c("2i_1", "2i_2", "ser_1", "ser_2") & nb_pmt_joshi_1 == 0)]
DT <- DT[!(condition_2 %in% c("2i_1", "2i_2", "ser_1", "ser_2") & nb_pmt_joshi_2 == 0)]
DT <- DT[!(condition_1 %in% c("sahlen") & nb_pmt_sahlen_1 == 0)]
DT <- DT[!(condition_2 %in% c("sahlen") & nb_pmt_sahlen_2 == 0)]

DT[, condition_1:=factor(condition_1, levels = c("2i_1", "2i_2", "ser_1", "ser_2", "sahlen"))]
DT[, condition_2:=factor(condition_2, levels = c("2i_1", "2i_2", "ser_1", "ser_2", "sahlen"))]

# Get spearman corr
r <- DT[, cor(reads_1, reads_2, method='spearman'), by=.(condition_1, condition_2)]$V1

# Remove diagonal entries
DT <- DT[condition_1 != condition_2]

# Make hexbin scatterplot
ggplot(DT, aes(x=reads_1, y=reads_2)) +
geom_hex(stat = "binhex", bins = 50) +
annotate(geom="text", x=3.5, y=9, label = sprintf("rho==%2.2f",round(r,2)), parse=T) +
facet_grid(condition_1 ~ condition_2) +
scale_x_continuous(breaks = c(0,5,10)) + 
scale_y_continuous(breaks = c(0,5,10)) + 
scale_fill_viridis() +
xlab("CHi-C library 1") +
ylab("CHi-C library 2") +
theme_SE() +
theme(legend.key.height = unit(15, "mm"))

# Save PDF
ggsave("IMG/FigS3B.pdf", width = 7, height = 5.5)
