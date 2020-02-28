library("GenomicRanges")

feature.matrix <- as.matrix(read.table(file="data/AllChIPseqFeatureMatrix.txt",sep="\t",head=T,row.names=1))

library(NMF)
png(file="ChIPseq_SElets_Heatmap_ESC_S7D.png",width=12,height=6,units='in',res=300)
aheatmap(feature.matrix[1:3,c(20,18,19)],Rowv=NA,Colv=NA,scale="col")
dev.off()

png(file="ChIPseq_SElets_Heatmap_ESC_S7E.png",width=12,height=6,units='in',res=300)
aheatmap(feature.matrix[1:3,c(16,17,15)],Rowv=NA,Colv=NA,scale="col")
dev.off()

png(file="ChIPseq_SElets_Heatmap_ESC_Fig5e.png",width=12,height=6,units='in',res=300)
# row order?
aheatmap(feature.matrix[1:3,c(20,3,6,7,8,4,0,1,2,9,10,11,12,14,21)+1],Rowv=NA,Colv=NA,scale="col")
dev.off()

# now Epi too

feature.matrixEpi <- as.matrix(read.table("data/AllChIPseqFeatureMatrix_Epi.txt",sep="\t",head=T,row.names=1))

png(file="ChIPseq_SElets_Heatmap_EpiLCplusA.png",width=8,height=4,units='in',res=300)
aheatmap(feature.matrixEpi[1:3,5:8],Rowv=NA,Colv=NA,scale="col")
dev.off()

png(file="ChIPseq_SElets_Heatmap_EpiSC.png",width=8,height=4,units='in',res=300)
aheatmap(feature.matrixEpi[1:3,c(9:12,17,16)],Rowv=NA,Colv=NA,scale="col")
dev.off()

# now add ChIP-exo analysis

chipexo.Rawcounts <- data.frame(chr=as.character(seqnames(enhancelet.gr)),start=as.numeric(start(enhancelet.gr)),end=as.numeric(end(enhancelet.gr)),class=featureType,Esrrb.GSE97304=rep(NA,length(enhancelet.gr)),Sox2.GSE54103=rep(NA,length(enhancelet.gr)),Stat3.GSE97304=rep(NA,length(enhancelet.gr)))

chipexo.rpkm <- data.frame(chr=as.character(seqnames(enhancelet.gr)),start=as.numeric(start(enhancelet.gr)),end=as.numeric(end(enhancelet.gr)),class=featureType,Esrrb.GSE97304=rep(NA,length(enhancelet.gr)),Sox2.GSE54103=rep(NA,length(enhancelet.gr)),Stat3.GSE97304=rep(NA,length(enhancelet.gr)))

chipexofiles <- c("GSE97304/bam/SRR5401052.bam","ChIPexo/bam/SRR1118256.bam","ChIPexo/bam/SRR5401053.bam")

for(i in 1:length(chipexofiles)){
	this.ChIP.reads <- readGAlignments(file=chipexofiles[i],param=ScanBamParam())
        cat("counting overlaps \n")
        this.ChIP.counts <- summarizeOverlaps(features=enhancelet.gr,reads=this.ChIP.reads)
	chipexo.Rawcounts[,i+4] <- assay(this.ChIP.counts)
	chipexo.rpkm[,i+4] <- 1e6*assay(this.ChIP.counts)/length(this.ChIP.reads)
}

chipexo.matrix <- t(sapply(c("PU","DM","PM"),function(x)colMeans(chipexo.rpkm[which(chipexo.rpkm$class==x),5:7])))
chipexo.matrix <- as.matrix(read.table("data/AllChIPseqFeatureMatrix_chipexo.txt",sep="\t",head=T,row.names=1))
			   
png(file="ChIPexo_SElets_Heatmap.png",width=6,height=6,units='in',res=300)
aheatmap(chipexo.matrix,Rowv=NA,Colv=NA,scale="col")
dev.off()
