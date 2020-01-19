######################################################################
#
# Plot 4C-seq interactions between Klf4 promoter and Klf4 SElets
#
# Wout Megchelenbrink
# Jan 17, 2020
######################################################################

rm(list=ls())
source("include/4cseq.inc.R")
source("include/style.R")

# Parameters
nb.dist         <- 1e6  # Near bait distance
normalization   <- "cis" # normalize using all reads in CIS
col.dark        <- c("green", "red", "blue")
min.rpkm        <- 1000 # Bins with less than 1000 RPKM are treated as background ligation and not considered for DE analysis
ymin            <- 0    # Lower bound of the y-axis  (x 1000 FPKM)
ymax            <- 13   # Upper bound of the y-axis (x 1000 FPKM)
XRANG           <- c(-90, 15) # Viewing window, 90kb upstream to 15kb downstream of the viewpoint

# Contrast for DEseq
contrasts <- list("2iL_vs_SL"=c("condition", "WT_2iL", "WT_SL"),
                  "KO_vs_SL"=c("condition", "KO_SL", "WT_SL"))

# Convert to 4Cseq format
fc <- fourcseq.to.4cseq(readRDS("DATA/4Cseq_SElets_KLF4.Rds"))

# Make figure that shows the interactions strength (running median of 7 neighboring Dpn2 frags)
run.med <- get.running.medians(fc,normalization = normalization, nb.dist = nb.dist, k=7, min.rpkm = min.rpkm)

(gg1 <- plot.4cseq(run.med, fc, min.rpkm = min.rpkm, 
                   xrange=XRANG, yrange=c(ymin, ymax), 
                   ybreaks = c(ymin, ymax), conditions = c("WT_2iL", "WT_SL", "KO_SL")))

tg <- tile.4cseq(fc, normalization = normalization, resolution = 5000, nb.dist = nb.dist)
DTX <- deseq2.4cseq(tg, run.med, contrasts = contrasts)

### DE interactions compared to SL-WT

# 2iL
DTX[above.threshold == T & abs(dist) < 1e5 & padj_2iL_vs_SL  < 0.05][order(padj_2iL_vs_SL)]

# SL KO
DTX[above.threshold == T & abs(dist) < 1e5 & padj_KO_vs_SL < 0.05][order(padj_KO_vs_SL)]

# Make fold change figure
(gg2 <- plot.4cseq.fc(DT.deseq = DTX, xrange=XRANG, nb.dist=nb.dist))

# Attach running medians and FC together and Save the PDF
export.4cseq.fig(gg1, gg2, vp, file.name = "IMG/Fig6B_KLF4_4Cseq.pdf")