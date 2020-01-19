######################################################################
# Wout Megchelenbrink
# Jan 18, 2020
#
# 4-way venn diagram of the overlap between the gene-SE interactions
# in Joshi, Sahlen, Novo (ChiC) and Dowen (ChiA-PET)
######################################################################

source("include/style.R")

DT <- fread("DATA/PMT_SE_interactions_quadruple_venn.tsv")
pdf("IMG/FigS3C.pdf", width = 4, height = 4)
VennDiagram::draw.quad.venn(area1 = nrow(DT[in.joshi == 1]),
                            area2 = nrow(DT[in.sahlen == 1]),
                            area3 = nrow(DT[in.cnovo == 1]),
                            area4 = nrow(DT[in.chiapet == 1]),
                            n12 =  nrow(DT[in.joshi == 1 & in.sahlen ==1]),
                            n13 = nrow(DT[in.joshi == 1 & in.cnovo ==1]),
                            n14 = nrow(DT[in.joshi == 1 & in.chiapet ==1]),
                            n23 = nrow(DT[in.sahlen == 1 & in.cnovo ==1]),
                            n24 = nrow(DT[in.sahlen == 1 & in.chiapet ==1]),
                            n34 = nrow(DT[in.cnovo == 1 & in.chiapet ==1]),
                            n123 = nrow(DT[in.joshi == 1 & in.sahlen ==1 & in.cnovo ==1]),
                            n124 = nrow(DT[in.joshi == 1 & in.sahlen ==1 & in.chiapet ==1]),
                            n234 = nrow(DT[in.sahlen == 1 & in.cnovo ==1 & in.chiapet ==1]),
                            n134 = nrow(DT[in.joshi == 1 & in.cnovo ==1 & in.chiapet ==1]),
                            n1234 = nrow(DT[in.joshi == 1 & in.sahlen ==1 & in.cnovo ==1 & in.chiapet ==1]),
                            fill = bs.col.dark[c(2,1,3,4)],
                            category = c("Joshi", "Sahlen", "Novo", "Dowen"))
dev.off()
# 30 + 36 + 52 + 108 = 226 shared between sahlen & joshi
