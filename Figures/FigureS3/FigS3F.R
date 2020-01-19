####################################################################
#
# Wout Megchelenbrink
# Jan. 17, 2020
# SE engaged in closest or more distal promoter interactions
###################################################################
rm(list=ls())

source("include/style.R")
chic <- unique(fread("DATA/CHIC_promoter_SE_interactions_Joshi_Sahlen_629_with_RPKM_ser_above_1_and_gene_status_KNOWN.tab")[gene_status == "KNOWN"], by=c("se_id","gene_name"))
se <- fread("DATA/Superenhancers.tsv")
DT <- merge(chic, se, by="se_id")


se.with.closest <- DT[gene_name == closest_expressed_gene & rpkm_ser >= 1, unique(se_id)]
se.with.distal  <- DT[gene_name != closest_expressed_gene & rpkm_ser >= 1, unique(se_id)]
se.with.both    <- intersect(se.with.closest, se.with.distal)
se.not.engaged  <- se[!se_id %in% chic$se_id,se_id]
se.not.expressed <- setdiff(chic[is.na(rpkm_ser) | rpkm_ser < 1, unique(se_id)], chic[rpkm_ser >= 1, unique(se_id)]) 

se[!se_id %in%  se.with.closest & !se_id %in% se.with.distal & !se_id %in% se.not.engaged & !se_id %in%  se.not.expressed]
length(se.with.closest) + length(se.with.distal) - length(se.with.both) + length(se.not.engaged) + length(se.not.expressed)

## Barplot of interaction categories
DT <- data.table(category=c("Closest and distal", "Only closest", "Only distal", "Only non-expressed", "No interaction"), 
                 N=c(length(se.with.both), length(se.with.closest)-length(se.with.both), length(se.with.distal)-length(se.with.both), length(se.not.engaged), length(se.not.expressed)))
DT[, category:=factor(category, levels = rev(c("Closest and distal", "Only closest", "Only distal", "Only non-expressed", "No interaction")))]

ggplot(DT, aes(x=category, y=N, fill=category, label=N)) +
geom_bar(stat = "identity") +
coord_flip() + 
scale_fill_manual(values = bs.col.dark[5:1], guide="none") +
theme_SE() +
xlab("") +
ylab("Interacting SE (#)") +
geom_text(nudge_x = 0, nudge_y = 10)

ggsave("IMG/FigS3F.pdf", width = 4.5, height = 2.5)


