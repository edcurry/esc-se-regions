######################################################################
#
#
# Wout Megchelenbrink
# January 18, 2020
#
# Sum 4C-seq interactions between each Klf4 associated SElet 
# and the KLF4 promoter in Esrrb WT and KO conditions
######################################################################

rm(list=ls())

source("include/4cseq.inc.R")
source("include/style.R")

# Params
nb.dist       <- 1e5    # Consider interactions within 100kb from the viewpoint (this is more than far enough for the Klf4 associated SE)
normalization <- "cis"  # normalize using all reads in CIS

# Contrast for DEseq
contrasts <- list("2iL_vs_SL"=c("condition", "WT_2iL", "WT_SL"),
                  "KO_vs_SL"=c("condition", "KO_SL", "WT_SL"))

fourcseq <- readRDS("DATA/4Cseq_SElets_KLF4.Rds")

# Convert
fc <- fourcseq.to.4cseq(fourcseq)

DT <- cbind(as.data.table(rowData(fc)), as.data.table(assay(fc, "counts")))
vp <- metadata(fc)$vp

FCseq <- normalize.4cseq(fc, normalization, nb.dist)
FCseq.gr <- makeGRangesFromDataFrame(FCseq)

# Get the KLF4 SElets and the control region
SElets <- fread("DATA/SElets.tsv")[chr == "chr4" & start > 55460000 & end < 55500000][type %in% c("INT", "DM", "PU"), .(SElet_id, type, selet_chr=chr, selet_start=start, selet_end=end)]
cntrl <- data.table(SElet_id=999, type="CTRL", selet_chr="chr4", selet_start=55503244-1000L, selet_end=55503244+2000L)
SElets <- rbind(SElets, cntrl)

# Convert to GRanges
SElets.gr <- makeGRangesFromDataFrame(SElets)

# Find overlap with 4C-seq fragments
ovl <- findOverlaps(SElets.gr, FCseq.gr)
DTX <- cbind(SElets[queryHits(ovl), ], FCseq[subjectHits(ovl), ])

# Sum interacting fragments per SElet
DTX <- DTX[, lapply(.SD, sum), by=.(SElet_id, type, selet_chr, selet_start, selet_end), .SDcols = c("WT_SL_1", "WT_SL_2", "WT_SL_3",
                                                                                                    "KO_SL_1", "KO_SL_2", "KO_SL_3")]

# To long format 
DTY <- melt.data.table(DTX, id.vars = c("SElet_id", "type"), measure.vars = c("WT_SL_1", "WT_SL_2", "WT_SL_3", "KO_SL_1", "KO_SL_2", "KO_SL_3"), variable.name = "genotype", value.name = "FPKM")
DTY[, genotype:=str_sub(genotype, end = 2)]

# Get mean and SEM
DTQ <- DTY[, .(FPKM.avg=mean(FPKM), FPKM.sem=sd(FPKM)/sqrt(3)), by=c("SElet_id", "type", "genotype")]
DTQ[, genotype:=factor(genotype, levels = c("WT", "KO"))]

# Make the plot
ggplot(DTQ, aes(x=paste0(SElet_id, "_", type), y=FPKM.avg, fill = genotype)) +
geom_bar(stat="identity", position = "dodge") +
geom_errorbar(aes(ymin=FPKM.avg-FPKM.sem, ymax=FPKM.avg+FPKM.sem), position = position_dodge2(width = .25, padding = .5)) +
xlab("SElet") +
ylab("4Cseq interaction (FPKM)") +
scale_fill_manual(values = c("#5082BE","#A6A6A6"))  +
theme_SE()

# Save PDF
ggsave("IMG/Fig6B.pdf", width = 5, height = 4)
