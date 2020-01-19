########################################################################
#
# Wout Megchelenbrink
# August 15, 2019
# Some (not to flexible and hard-coded) helper functions for the 4Cseq
#
########################################################################

source("include/style.R")


fourcseq.to.4cseq.list <- function(lst)
{
  lst.out <- list()
  
  for(i in 1:length(lst))
  {
    lst.out[[i]] <- fourcseq.to.4cseq(lst[[i]])  
  }
  
  return(lst.out)
}


fourcseq.to.4cseq <- function(fourcseq)
{
  cnts <- counts(fourcseq)
  vp <- as.character(colData(fourcseq)$viewpoint[1])
  colnames(cnts) <- str_replace(colnames(cnts), paste0(vp, "_"), "")
  fcs <- SummarizedExperiment(assays = list("counts" = cnts))
  fcs@metadata$vp <- vp
  fcs@metadata$vp.chr <- as.character(colData(fourcseq)$chr[1])
  fcs@metadata$vp.coord <- as.integer(colData(fourcseq)$primerStart[1] + colData(fourcseq)$primerEnd[1])/2
  rows <- as.data.table(rowRanges(fourcseq))[,1:3]
  setnames(rows,  c("chr","start", "end"))
  rowData(fcs) <- rows
  
  cd <- as.data.table(colData(fourcseq)[ ,c("condition", "replicate")])
  cd[, replicate:=factor(replicate)]
  cd[, condition:=factor(condition)]
  
  colData(fcs) <- DataFrame(row.names = colnames(cnts), cd)
  return(fcs)
}


get.running.medians <- function(fc, normalization = "wg", nb.dist = 1e6, k=7, min.rpkm = 250)
{
  # Get the sample names
  cnt <- assay(fc, "counts")
  samples <- colnames(cnt)
  
  # Bind the rowranges tot the count data and compute distance to viewpoint
  fcseq <- cbind(as.data.table(rowData(fc)), as.data.table(cnt))
  fcseq[, vp.chr:=metadata(fc)$vp.chr]
  fcseq[, vp.coord:=metadata(fc)$vp.coord]
  fcseq[, dist:=as.integer((start+end)/2 - vp.coord)]
  fcseq[, width:=end-start]

  if(normalization == "cis")
  {
    cat("doing cis normalization\n")
    fcseq <- fcseq[chr == vp.chr]
    
  } else if(normalization == "nb") # near bait
  {
    cat("doing NB normalization\n")
    fcseq <- fcseq[chr == vp.chr & abs(dist) <= nb.dist]
  }
  
  # Normalize reads
  cat("now normalizing with DESeq2\n")
  
  cnt <- fcseq[, -c("chr", "start", "end", "vp.chr", "vp.coord", "dist", "width"), with=F]
  
  dds <- DESeqDataSetFromMatrix(countData = cnt, colData = colData(fc), design = ~condition)
  dds <- estimateSizeFactors(dds)
  #dds <- DESeq(dds, parallel = F, fitType = "local", test = "Wald") # Fittype = local makes quite a difference
  
  mcols(dds)$basepairs <- fcseq[, end-start]
  rpkm.scaled <- DESeq2::fpkm(dds, robust = T)
  fcseq <- cbind(fcseq[, .(chr, start, end, vp.chr, vp.coord, dist)], rpkm.scaled)
  
  # Get the near bait 4c-seq
  fcseq.smooth <- apply(data.matrix(fcseq[, samples, with=F]), 2, function(x) rollmean(x, k, "extend"))
  
  #### Annotate which fragments have signal above the threshold
  libs.WT.2iL <- samples[grep("WT_2iL", samples)]
  libs.WT.SL <- samples[grep("WT_SL", samples)]
  libs.KO.SL <- samples[grep("KO_SL", samples)]
  
  # Normalized counts
  ID <- which(matrixStats::rowMedians(fcseq.smooth[, libs.WT.2iL]) >= min.rpkm | 
                matrixStats::rowMedians(fcseq.smooth[, libs.WT.SL]) >= min.rpkm | 
                matrixStats::rowMedians(fcseq.smooth[, libs.KO.SL]) >= min.rpkm)
  
  fcseq.smooth <- cbind(fcseq[, .(chr, start, end, vp.chr, vp.coord, dist)], as.data.table(fcseq.smooth))
  fcseq.smooth[, above.threshold:=F]
  fcseq.smooth[ID, above.threshold:=T]
  
  return(fcseq.smooth)
}


# Bin the genome into tiles of X bp
getBinnedGenome <- function(resolution = 5000, chr = NULL)
{
  sls <- seqinfo(BSgenome.Mmusculus.UCSC.mm9)@seqlengths
  names(sls) <- seqinfo(BSgenome.Mmusculus.UCSC.mm9)@seqnames
  
  if(!is.null(chr))
  {
    # Tile width cannot exceed chromosome length
    resolution <- pmin(sls[[chr]], resolution)
    
    cat(sprintf("BINNING CHR :: %s\n", chr))
    seq <- sls[[chr]]
    names(seq) <- chr
    binned.genome <- tileGenome(seq, tilewidth=resolution, cut.last.tile.in.chrom=TRUE) 
    print(length(binned.genome))
  } else
  {
    binned.genome <- tileGenome(sls, tilewidth=resolution, cut.last.tile.in.chrom=TRUE) 
  }
  
  binned.genome <- as.data.table(binned.genome)[, 1:3]
  setnames(binned.genome, c("chr", "start", "end"))
  
  return(binned.genome)
}


normalize.4cseq <- function(fc, normalization = "wg", nb.dist = 1e6)
{
  # Get the sample names
  cnt <- assay(fc, "counts")
  samples <- colnames(cnt)
  
  # Bind the rowranges tot the count data and compute distance to viewpoint
  fcseq <- cbind(as.data.table(rowData(fc)), as.data.table(cnt))
  fcseq[, vp.chr:=metadata(fc)$vp.chr]
  fcseq[, vp.coord:=metadata(fc)$vp.coord]
  fcseq[, dist:=as.integer((start+end)/2 - vp.coord)]
  fcseq[, width:=end-start]
  
  if(normalization == "cis")
  {
    cat("doing cis normalization\n")
    fcseq <- fcseq[chr == vp.chr]
    
  } else if(normalization == "nb") # near bait
  {
    cat("doing NB normalization\n")
    fcseq <- fcseq[chr == vp.chr & abs(dist) <= nb.dist]
  }
  
  # Normalize reads
  cat("now normalizing with DESeq2\n")
  
  cnt <- fcseq[, -c("chr", "start", "end", "vp.chr", "vp.coord", "dist", "width"), with=F]
  
  dds <- DESeqDataSetFromMatrix(countData = cnt, colData = colData(fc), design = ~condition)
  dds <- estimateSizeFactors(dds)
  #dds <- DESeq(dds, parallel = F, fitType = "local", test = "Wald") # Fittype = local makes quite a difference
  
  mcols(dds)$basepairs <- fcseq[, end-start]
  rpkm.scaled <- DESeq2::fpkm(dds, robust = T)
  fcseq <- cbind(fcseq[, .(chr, start, end, vp.chr, vp.coord, dist)], rpkm.scaled)

  return(fcseq)
}


# Bin the genome into tiles of X bp
getBinnedGenome <- function(resolution = 5000, chr = NULL)
{
  sls <- seqinfo(BSgenome.Mmusculus.UCSC.mm9)@seqlengths
  names(sls) <- seqinfo(BSgenome.Mmusculus.UCSC.mm9)@seqnames
  
  if(!is.null(chr))
  {
    # Tile width cannot exceed chromosome length
    resolution <- pmin(sls[[chr]], resolution)
    
    cat(sprintf("BINNING CHR :: %s\n", chr))
    seq <- sls[[chr]]
    names(seq) <- chr
    binned.genome <- tileGenome(seq, tilewidth=resolution, cut.last.tile.in.chrom=TRUE) 
    print(length(binned.genome))
  } else
  {
    binned.genome <- tileGenome(sls, tilewidth=resolution, cut.last.tile.in.chrom=TRUE) 
  }
  
  binned.genome <- as.data.table(binned.genome)[, 1:3]
  setnames(binned.genome, c("chr", "start", "end"))
  
  return(binned.genome)
}






tile.4cseq <- function(fcs, normalization, resolution = 5000, nb.dist=1e6)
{
  cnts <- cbind(as.data.table(rowData(fcs)), as.data.table(assay(fcs, "counts")))

  # Bin the genome
  bg <- getBinnedGenome(resolution)
  #bg[chr ==  metadata(fcs)$vp.chr]
 
  bg[, vp.chr:=metadata(fcs)$vp.chr]
  bg[, vp.coord:=metadata(fcs)$vp.coord]
  
   bg[chr==vp.chr, dist:=as.integer((start+end)/2 - vp.coord)]
  #  ID <- bg[chr==vp.chr, .I[abs(dist) <= min(abs(dist))]]
  # 
  # foo <- bg[chr==vp.chr][ID, dist]
  # 
  # if(length(foo) != 1)
  # {
  #   print(foo)
  #   stop("Length foo != 1")
  # }
# 
#   type(foo)
#   bg[chr==vp.chr,]
#   bg[order(abs(dist))]
#   bg[chr==vp.chr, start:=start-foo]
#   bg[chr==vp.chr, end:=end-foo]
# 
#   bg[chr==vp.chr, dist:=as.integer((start+end)/2 - vp.coord)]
#   bg[chr == vp.chr, min(abs(dist))]

  
    #  bg <- getBinnedGenome(resolution, chr = metadata(fcs)$vp.chr)
  setnames(bg, c("bin_chr", "bin_start", "bin_end", "vp.chr", "vp.coord", "dist"))

  # Make GRanges of counts and bins
  bg.gr <- makeGRangesFromDataFrame(bg, seqnames.field = "bin_chr", start.field = "bin_start", end.field = "bin_end")
  counts.gr <- makeGRangesFromDataFrame(cnts)
  
  # Find overlaps and concatenate
  ovl <- findOverlaps(counts.gr, bg.gr)
  overlaps <- pintersect(counts.gr[queryHits(ovl)], bg.gr[subjectHits(ovl)])
  pct.ovl <- width(overlaps) / width(bg.gr[subjectHits(ovl)]) * 100
  DT <- cbind(cnts[queryHits(ovl), ], bg[subjectHits(ovl),], pct.ovl)
  
  # Assign to bin with most overlap
  # DT <- DT[DT[, .I[pct.ovl >= max(pct.ovl)], by=.(chr, start, end)]$V1,]
  
  DT[chr == vp.chr, min(abs(dist))]
  # Sum reads per tile
  DTX <- DT[, lapply(.SD, sum), by=.(bin_chr, bin_start, bin_end, vp.chr, vp.coord, dist), .SDcols=colnames(fcs)]
 # DTX[, vp.chr:=metadata(fcs)$vp.chr]
#  DTX[, vp.coord:=metadata(fcs)$vp.coord]
#  DTX[, dist:=(bin_start+bin_end)/2 - vp.coord]

  # change this  .. 
  
  if(normalization == "cis")
  {
    cat("selecting CIS tiles for DESeq2\n")
    DTX <- DTX[bin_chr == vp.chr]
  } else if(normalization == "nb") # near bait
  {
    cat("selecting NB tiles for DEseq2\n")
    DTX <- DTX[bin_chr == vp.chr & abs(dist) <= nb.dist]
  }
  
  # Create new SE object
  fcs.out <- SummarizedExperiment(assays = list("counts"=data.matrix(DTX[, colnames(fcs), with=F])))
  rowData(fcs.out) <- DTX[, .(vp.chr, vp.coord, chr=bin_chr, start=bin_start, end=bin_end, dist)]
  metadata(fcs.out) <- metadata(fcs)
  colData(fcs.out) <- colData(fcs)
  
  return(fcs.out)
}

deseq2.4cseq <- function(tg, run.med, contrasts = list("ser_vs_2i"=c("condition", "ser", "2i")))
{
  cnts <- round(assay(tg, "counts"))
  dds <- DESeqDataSetFromMatrix(countData = cnts, colData = colData(tg), design = ~condition)
  dds <- estimateSizeFactors(dds)
  dds <- DESeq(dds, parallel = F, fitType = "local", test = "Wald") # Fittype = local makes quite a difference
  
  mcols(dds)$basepairs <- rowData(tg)[, "end"]-rowData(tg)[, "start"]
  rpkm <- DESeq2::fpkm(dds, robust = T)
  
  DTX <- cbind(as.data.table(rowData(tg)), rpkm)

  lst <- list()
  for(i in 1:length(contrasts))
  {
    contrast <- contrasts[[i]]
    res <- as.data.table(results(dds, contrast=contrast))
    res <- res[, .(baseMean, log2FoldChange, lfcSE, pvalue, padj)]
    setnames(res, paste0(c("baseMean_","log2fc_", "lfcSE_", "pval_", "padj_"), names(contrasts)[[i]]))
    DTX <- cbind(DTX, res)    
  }
  
  runmed.gr <- makeGRangesFromDataFrame(run.med[above.threshold == T])
  DTX.gr <- makeGRangesFromDataFrame(DTX)
  
  ovl <- findOverlaps(DTX.gr, runmed.gr)
  DTX[, above.threshold:=F]
  DTX[unique(queryHits(ovl)), above.threshold:=T]
  
  return(DTX)
}

plot.4cseq <- function(run.med, fc, min.rpkm = 250, xrange=c(-100, 100), yrange=c(0,5.5), ybreaks=c(0,5), 
                       conditions = c("WT_2iL", "WT_SL", "KO_SL"))
{
  run.med <- run.med[dist >= xrange[1]*1000 & dist <= xrange[2]*1000]
  DTX <- melt.data.table(run.med[, c("vp.chr", "vp.coord", "chr", "start", "end", "dist", colnames(assay(fc, "counts"))), with=F],
                         id.vars = c("chr", "start", "end", "vp.chr", "vp.coord", "dist"))
  DTX[grep("WT_2iL", variable), condition:="WT_2iL"]
  DTX[grep("WT_SL", variable), condition:="WT_SL"]
  DTX[grep("KO_SL", variable), condition:="KO_SL"]
  
  DTX[, condition:=factor(condition, levels=conditions)]
  
  DTX <- DTX[, .(avg=median(value), stdev=sd(value)), by=.(chr, start, end, dist, condition)]
  DTX[, sem:=stdev/sqrt(3)]
  DTX[, avg:=avg/1000]
  DTX[, sem:=sem/1000]
  DTX[, stdev:=stdev/1000]
  DTX[, dist:=dist/1000]
 
   # Get viewing window
  
  # Why do I do this??
  #DTX[avg < min.rpkm/1000, avg:=-1000]

  # Make plot
   gg <- 
     ggplot(DTX, aes(x=dist, y=avg, color=condition, fill=condition)) +
        geom_line() + #(stat = 'identity') +
      geom_ribbon(aes(x=dist, ymin=avg - sem, ymax=avg + sem), color = NA) +
        geom_hline(yintercept = min.rpkm/1000, linetype = "dashed", color = "#222222", cex = .2) +
        scale_x_continuous(labels = comma, expand = c(0, 0)) +
        scale_y_continuous(limits = yrange, breaks = ybreaks, labels = comma,  expand = c(0, 0)) +
        scale_color_manual(values = colors.dark[2:4]) +
         scale_fill_manual(values = alpha(colors.dark[2:4], .35)) +
             ylab(NULL) +
              theme_SE() +
             guides(fill="none", color="none") +
             facet_wrap(~condition, nrow = DTX[, length(unique(condition))], strip.position="right") +
             theme(strip.text = element_text(size = 5),
                   strip.switch.pad.grid = unit(0.5, "lines"),
                   strip.background = element_blank(),
                   panel.spacing = unit(0.5, "lines"),
                   axis.line.x = element_line(size = .75),
                   axis.line.y = element_line(size = .5),
                   axis.line.x.bottom = element_blank(),
                   axis.line.y.left = element_line(size = .25))

  return(gg)
}



plot.4cseq.fc <- function(DT.deseq, xrange, nb.dist)
{
  DT.deseq[, binID:=.I]
  DT.deseq <- DT.deseq[abs(dist) <= nb.dist]
  DT.deseq <- melt.data.table(DT.deseq, id.vars = c("binID", "vp.chr", "vp.coord", "chr", "start", "end", "dist", "above.threshold"),
                              measure.vars = list(c("log2fc_2iL_vs_SL", "log2fc_KO_vs_SL"),
                                                  c("lfcSE_2iL_vs_SL", "lfcSE_KO_vs_SL"),
                                                  c("padj_2iL_vs_SL", "padj_KO_vs_SL")),
                              variable.name = "condition", value.name = c("log2fc", "log2fcSE", "padj"))

  DT.deseq[, condition:=ifelse(condition == 1, "2iL_vs_SL", "KO_vs_SL")] 
  
  # bins <- DT.deseq[above.threshold == T & padj < 0.05, unique(binID)]
  # DT.deseq <- DT.deseq[!binID %in% bins, padj:=1]
  # DT.deseq <- DT.deseq[!binID %in% bins, log2fc:=0]
   DT.deseq[, score:=-log10(padj)]
   DT.deseq[, dist:=dist/1000]
  # 
  DT.deseq <- DT.deseq[dist > xrange[1] & dist <= xrange[2]]
  
  DT.deseq <-  DT.deseq[abs(dist) >= 10]
  
  # Significance marks
  DT.deseq[padj < .1, label:="."]
  DT.deseq[padj < .05, label:="*"]
  DT.deseq[padj < .001, label:="**"]
  DT.deseq[padj < 1e-5, label:="***"]
  
  DT.deseq[, condition:=factor(condition, levels=c("2iL_vs_SL", "KO_vs_SL"))]
  DT.deseq[above.threshold == F, log2fc:=NA]
  DT.deseq[above.threshold == F, label:=NA]
  DT.deseq[, fc.sign:=factor(sign(log2fc), levels = c(-1,0,1))]
  
  
  min.fc <- floor(DT.deseq[, min(log2fc, na.rm = T)])
  max.fc <- ceiling(DT.deseq[, max(log2fc,  na.rm = T)])
  #   
  #DT.deseq <- DT.deseq[above.threshold == T & padj < .1]
  
  gg <- ggplot(DT.deseq, aes(x=dist, y=log2fc,  color = condition, label=label)) +
    #geom_bar(stat = "identity", position = position_dodge(width = .9)) +
    geom_errorbar(aes(ymin = log2fc - log2fcSE, ymax=log2fc + log2fcSE), size=.15) +
    geom_point(cex=.1) +
    scale_color_manual(values = c(bs.col.dark[3], bs.col.dark[4])) +
   # scale_fill_manual(values = c(bs.col.light[1], "purple", bs.col.light[2])) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "black", cex = .1, linetype="dotted") +
    xlab("Genomic distance (kb)") +
    ylab(NULL) +
    geom_text(data = DT.deseq[fc.sign==1], nudge_y=0, color = "black", cex=4) +
    geom_text(data = DT.deseq[fc.sign==-1], nudge_y=0, color = "black",  cex=4) +
    scale_x_continuous(expand = c(0, 0), limits = c(xrange[1], xrange[2])) +
    scale_y_continuous(expand = c(0, 0), limits = c(min.fc -.5, max.fc + .5), breaks = c(min.fc, max.fc)) +
    theme_SE() +
    guides(color="none", fill = "none") + 
    facet_wrap(~condition, nrow = 2, strip.position="right") +
    theme(strip.text = element_text(size = 5),
          strip.switch.pad.grid = unit(0.5, "lines"),
          strip.background = element_blank(), 
           panel.spacing = unit(0.5, "lines"),
           axis.line.x = element_line(size = .75),
           axis.line.y = element_line(size = .5),
           axis.line.x.bottom = element_line(size = .25),
           axis.line.y.left = element_line(size = .25))

  return(gg)
}





export.4cseq.fig <- function(gg1, gg2, vp, file.name)
{
  gg1 <- gg1 + theme(axis.text.x = element_blank(),
                     axis.title.x = element_blank(), 
                     axis.ticks.x = element_blank())
  
  gg1 <- ggplotGrob(gg1)
  gg2 <- ggplotGrob(gg2)
  maxwidths <- grid::unit.pmax(gg1$widths[2:5], gg2$widths[2:5])
  gg1$widths[2:5] <- as.list(maxwidths)
  gg2$widths[2:5] <- as.list(maxwidths)
  
  blank <- rectGrob(gp=gpar(col="white")) # make a white spacer grob
  g <- gtable_matrix(name = "demo",
                     grobs = matrix(list(gg1, blank, gg2), nrow = 3), 
                     widths = unit(10, "cm"),
                     heights = unit(c(4, .1, 4), "cm"))

  grid.newpage()
  pdf(file.name, width = 10, height = 4)
  grid.draw(g)
  dev.off()
}
