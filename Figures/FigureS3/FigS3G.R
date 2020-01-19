#################################################################
#
# Wout Megchelenbrink
# Jan 17, 2020
#
# CHiC read intensity at closest and distal SE connected genes.
#################################################################

rm(list=ls())
source("include/style.R")

chic <- unique(fread("DATA/CHIC_promoter_SE_interactions_Joshi_Sahlen_629_with_RPKM_ser_above_1_and_gene_status_KNOWN.tab"), by=c("se_id","gene_name"))[gene_status == "KNOWN" & rpkm_ser >= 1]
se <- fread("DATA/Superenhancers.tsv")

DT <- merge(chic, se, by="se_id")

DTX <- rbind(DT[reads_sahlen >= 5 & score_sahlen >= 5,  .(se_id, gene_name, type="Sahlen", is.closest = ifelse(gene_name == closest_expressed_gene, "Closest expressed", "Distal"), chic_reads = reads_sahlen)],
             DT[sqrt(reads_ser_1*reads_ser_2) >= 5 & score_ser >= 5, .(se_id, gene_name, type="Joshi",is.closest = ifelse(gene_name == closest_expressed_gene, "Closest expressed", "Distal"),  chic_reads = sqrt(reads_ser_1*reads_ser_2))])

ggplot(DTX, aes(x=type, y=log2(chic_reads+1), fill = is.closest)) +
geom_boxplot(notch = T) +
scale_fill_manual(values = bs.col.dark[2:1]) + 
xlab("") +
ylab("C-HiC interaction strength (log2 reads)") +
theme_SE() +
theme(legend.position = "bottom", legend.title = element_blank())

ggsave(filename = "IMG/FigS3G.pdf", width = 3.5, height = 4.5)



