################################################################
# 
# Wout Megchelenbrink
# Jan 17, 2020
#
# Promoter-SE interactions
###############################################################

source("include/style.R")
chic <- fread("DATA/CHIC_promoter_SE_interactions_Joshi_Sahlen_629_with_RPKM_ser_above_1_and_gene_status_KNOWN.tab")[rpkm_ser >= 1 & gene_status == "KNOWN"]

chic <- chic[, lapply(.SD, max), by=.(gene_name, se_id), .SDcols=c("reads_ser", "score_ser", "reads_sahlen", "score_sahlen")]

# Take only interactions with expressed promoters
chic[reads_ser >= 5 & score_ser >= 5, in.joshi:=1]
chic[reads_sahlen >= 5 & score_sahlen >= 5, in.sahlen:=1]

chic[is.na(in.joshi), in.joshi:=0]
chic[is.na(in.sahlen), in.sahlen:=0]
chic[, .N, by=.(in.joshi, in.sahlen)]


pdf("IMG/FigS3D.pdf", width = 4, height = 4)
VennDiagram::draw.pairwise.venn(area1 = nrow(chic[in.joshi == 1]),
                                area2 = nrow(chic[in.sahlen == 1]),
                                cross.area = nrow(chic[in.joshi == 1 & in.sahlen == 1]),
                                fill = bs.col.dark[2:1],
                                category = c("J", "S"))
dev.off()