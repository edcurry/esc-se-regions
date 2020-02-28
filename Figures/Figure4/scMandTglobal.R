# load PM/DU/PU regions

library(GenomicRanges)
library(gplots)

PM.bed <- read.table("data/PM_mm10.bed",sep="\t",head=F)
PM.gr <- GRanges(seqnames=gsub(as.character(PM.bed[[1]]),pattern="chr",replace=""),ranges=IRanges(start=PM.bed[[2]],end=PM.bed[[3]]))

DM.bed <- read.table("data/DMR_mm10.bed",sep="\t",head=F)
DM.gr <- GRanges(seqnames=gsub(as.character(DM.bed[[1]]),pattern="chr",replace=""),ranges=IRanges(start=DM.bed[[2]],end=DM.bed[[3]]))

PU.bed <- read.table("data/PU_mm10.bed",sep="\t",head=F)
PU.gr <- GRanges(seqnames=gsub(as.character(PU.bed[[1]]),pattern="chr",replace=""),ranges=IRanges(start=PU.bed[[2]],end=PU.bed[[3]]))

## load each sample, for each region get M & U reads, store (M+U) & M/U
#
getMethRegion <- function(x,y){
	counts <- subsetByOverlaps(y,x)
	sum(counts$M)/(sum(counts$M)+sum(counts$U))
}

getCovRegion <- function(x,y){
        counts <- subsetByOverlaps(y,x)
        sum(counts$M+counts$U)
}

# obtain coverage files from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE68642&format=file
# for speed of processing, CpG-level read count files were subsetted to regions of interest using 'bedtools intersect', resulting in a directory of files corresponding to each region class (DMR, PU or PM)

all.dm.files <- list.files("data/DMRcov_GSE68642_mm10")
all.dm.gr <- list()
for(i in 1:length(all.dm.files)){
	cat(paste("getting DM methylation for cell",i,"of",length(all.dm.files),"\n"))
	this.bed <- read.table(paste("data/DMRcov_GSE68642_mm10/",all.dm.files[i],sep=""))
	all.dm.gr[[i]] <- subsetByOverlaps(GRanges(seqnames=as.character(this.bed[[1]]),ranges=IRanges(start=this.bed[[2]],end=this.bed[[3]]),M=this.bed[[5]],U=this.bed[[6]]),DM.gr)
}

all.pm.files <- list.files("data/PMcov_GSE68642_mm10")
all.pm.gr <- list()
for(i in 1:length(all.pm.files)){
        cat(paste("getting PM methylation for cell",i,"of",length(all.pm.files),"\n"))
        this.bed <- read.table(paste("data/PMcov_GSE68642_mm10/",all.pm.files[i],sep=""))
        all.pm.gr[[i]] <- subsetByOverlaps(GRanges(seqnames=as.character(this.bed[[1]]),ranges=IRanges(start=this.bed[[2]],end=this.bed[[3]]),M=this.bed[[5]],U=this.bed[[6]]),PM.gr)
}

all.pu.files <- list.files("PUcov_GSE68642_mm10")
all.pu.gr <- list()
for(i in 1:length(all.pu.files)){
        cat(paste("getting PU methylation for cell",i,"of",length(all.pu.files),"\n"))
        this.bed <- read.table(paste("data/PUcov_GSE68642_mm10/",all.pu.files[i],sep=""))
        all.pu.gr[[i]] <- subsetByOverlaps(GRanges(seqnames=as.character(this.bed[[1]]),ranges=IRanges(start=this.bed[[2]],end=this.bed[[3]]),M=this.bed[[5]],U=this.bed[[6]]),PU.gr)
}

# _overall_ analysis for each region type

full.countmatrix <- array(NA,dim=c(4,length(all.dm.gr)))
rownames(full.countmatrix) <- c("dm","pm","pu","cov")
colnames(full.countmatrix) <- sapply(strsplit(all.pu.files,split=".",fixed=T),function(x)x[[1]][1])
for(i in 1:length(all.dm.gr)){
full.countmatrix[1,i] <- sum(all.dm.gr[[i]]$M)/(sum(all.dm.gr[[i]]$M)+sum(all.dm.gr[[i]]$U))
full.countmatrix[2,i] <- sum(all.pm.gr[[i]]$M)/(sum(all.pm.gr[[i]]$M)+sum(all.pm.gr[[i]]$U))
full.countmatrix[3,i] <- sum(all.pu.gr[[i]]$M)/(sum(all.pu.gr[[i]]$M)+sum(all.pu.gr[[i]]$U))
full.countmatrix[4,i] <- sum(all.dm.gr[[i]]$M)+sum(all.dm.gr[[i]]$U)+sum(all.pm.gr[[i]]$M)+sum(all.pm.gr[[i]]$U)+sum(all.pu.gr[[i]]$M)+sum(all.pu.gr[[i]]$U)
}
full.countmatrix[4,] <- (full.countmatrix[4,]-min(full.countmatrix[4,]))/(max(full.countmatrix[4,])-min(full.countmatrix[4,]))

dm.cellCodeCols = rep("purple",length(all.dm.files))
dm.cellCodeCols[grep("2I",all.dm.files)] <- "blue"

png(file="EnhanceletRegionMeth_ByCell.png",height=900)
heatmap.2(t(full.countmatrix[1:3,]),scale="none",trace="none",RowSideColors=dm.cellCodeCols,col=bluered(100),Colv=NA)
dev.off()

dmr.clusters = cutree(hclust(dist(full.countmatrix[1,])),k=2)
