library("GenomicRanges")

chipexo.Rawcounts <- data.frame(chr=as.character(seqnames(enhancelet.gr)),start=as.numeric(start(enhancelet.gr)),end=as.numeric(end(enhancelet.gr)),class=featureType,Esrrb.GSE97304=rep(NA,length(enhancelet.gr)),Sox2.GSE54103=rep(NA,length(enhancelet.gr)),Stat3.GSE97304=rep(NA,length(enhancelet.gr)))

chipexo.rpkm <- data.frame(chr=as.character(seqnames(enhancelet.gr)),start=as.numeric(start(enhancelet.gr)),end=as.numeric(end(enhancelet.gr)),class=featureType,Esrrb.GSE97304=rep(NA,length(enhancelet.gr)),Sox2.GSE54103=rep(NA,length(enhancelet.gr)),Stat3.GSE97304=rep(NA,length(enhancelet.gr)))

chipexofiles <- c("data/SRR5401052.bam","data/SRR1118256.bam","data/SRR5401053.bam")

for(i in 1:length(chipexofiles)){
	this.ChIP.reads <- readGAlignments(file=chipexofiles[i],param=ScanBamParam())
        cat("counting overlaps \n")
        this.ChIP.counts <- summarizeOverlaps(features=enhancelet.gr,reads=this.ChIP.reads)
	chipexo.Rawcounts[,i+4] <- assay(this.ChIP.counts)
	chipexo.rpkm[,i+4] <- 1e6*assay(this.ChIP.counts)/length(this.ChIP.reads)
}

chipexo.matrix <- t(sapply(c("PU","DM","PM"),function(x)colMeans(chipexo.rpkm[which(chipexo.rpkm$class==x),5:7])))

library(NMF)
png(file="ChIPexo_SElets_Heatmap.png",width=6,height=6,units='in',res=300)
aheatmap(chipexo.matrix,Rowv=NA,Colv=NA,scale="col")
dev.off()
