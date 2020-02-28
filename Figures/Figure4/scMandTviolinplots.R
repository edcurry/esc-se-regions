# load PM/DU/PU regions

library(GenomicRanges)

PM.bed <- read.table("PM_mm10.bed",sep="\t",head=F)
PM.gr <- GRanges(seqnames=gsub(as.character(PM.bed[[1]]),pattern="chr",replace=""),ranges=IRanges(start=PM.bed[[2]],end=PM.bed[[3]]))

DM.bed <- read.table("DMR_mm10.bed",sep="\t",head=F)
DM.gr <- GRanges(seqnames=gsub(as.character(DM.bed[[1]]),pattern="chr",replace=""),ranges=IRanges(start=DM.bed[[2]],end=DM.bed[[3]]))

PU.bed <- read.table("PU_mm10.bed",sep="\t",head=F)
PU.gr <- GRanges(seqnames=gsub(as.character(PU.bed[[1]]),pattern="chr",replace=""),ranges=IRanges(start=PU.bed[[2]],end=PU.bed[[3]]))

## load each sample, for each region get M & U reads, store (M+U) & M/U
#
getMethRegion <- function(x,y){
	counts <- subsetByOverlaps(y,ranges=x)
	sum(counts$M)/(sum(counts$M)+sum(counts$U))
}

getCovRegion <- function(x,y){
        counts <- subsetByOverlaps(y,ranges=x)
        sum(counts$M+counts$U)
}

all.dm.files <- list.files("DMRcov_GSE68642_mm10")
all.dm.meth <- list()
all.dm.cov <- list()
for(i in 1:length(all.dm.files)){
	cat(paste("getting DM methylation for cell",i,"of",length(all.dm.files),"\n"))
	this.bed <- read.table(paste("DMRcov_GSE68642_mm10/",all.dm.files[i],sep=""))
	this.gr <- GRanges(seqnames=as.character(this.bed[[1]]),ranges=IRanges(start=this.bed[[2]],end=this.bed[[3]]),M=this.bed[[5]],U=this.bed[[6]])
	all.dm.meth[[i]] <- sapply(DM.gr,getMethRegion,y=this.gr)
	all.dm.cov[[i]] <- sapply(DM.gr,getCovRegion,y=this.gr)
}

all.pm.files <- list.files("PMcov_GSE68642_mm10")
all.pm.meth <- list()
all.pm.cov <- list()
for(i in 1:length(all.pm.files)){
        cat(paste("getting PM methylation for cell",i,"of",length(all.pm.files),"\n"))
        this.bed <- read.table(paste("PMcov_GSE68642_mm10/",all.pm.files[i],sep=""))
	this.gr <- GRanges(seqnames=as.character(this.bed[[1]]),ranges=IRanges(start=this.bed[[2]],end=this.bed[[3]]),M=this.bed[[5]],U=this.bed[[6]])
	all.pm.meth[[i]] <- sapply(PM.gr,getMethRegion,y=this.gr)
	all.pm.cov[[i]] <- sapply(PM.gr,getCovRegion,y=this.gr)
}

all.pu.files <- list.files("PUcov_GSE68642_mm10")
all.pu.meth <- list()
all.pu.cov <- list()
for(i in 1:length(all.pu.files)){
        cat(paste("getting PU methylation for cell",i,"of",length(all.pu.files),"\n"))
        this.bed <- read.table(paste("PUcov_GSE68642_mm10/",all.pu.files[i],sep=""))
	this.gr <- GRanges(seqnames=as.character(this.bed[[1]]),ranges=IRanges(start=this.bed[[2]],end=this.bed[[3]]),M=this.bed[[5]],U=this.bed[[6]])
        all.pu.meth[[i]] <- sapply(PU.gr,getMethRegion,y=this.gr)
	all.pu.cov[[i]] <- sapply(PU.gr,getCovRegion,y=this.gr)
}

dm.cellCodeCols = rep("purple",length(all.dm.files))
dm.cellCodeCols[grep("2I",all.dm.files)] <- "blue"
all.celltypes <- c('2I','serum')[2-as.numeric(dm.cellCodeCols=='blue')]

# produce violin plots for DM, PU (& PM?) for each cell
# => make dataframe with one row per region,cell. colnames: cell, region.type, meth
allcells <- sapply(strsplit(all.pu.files,split=".",fixed=T),function(x)x[[1]][1])
allcells.names <- c()
allcells.meth <- c()
allcells.regions <- c()
allcells.types <- c()
for(i in 1:length(allcells)){
	allcells.meth <- c(allcells.meth,all.dm.meth[[i]][which(all.dm.cov[[i]]>4)])
	allcells.regions <- c(allcells.regions,rep('DM',length(which(all.dm.cov[[i]]>4))))
	allcells.names <- c(allcells.names,rep(allcells[i],length(which(all.dm.cov[[i]]>4))))
	allcells.types <- c(allcells.types,rep(all.celltypes[i],length(which(all.dm.cov[[i]]>4))))
	allcells.meth <- c(allcells.meth,all.pm.meth[[i]][which(all.pm.cov[[i]]>4)])
	allcells.regions <- c(allcells.regions,rep('PM',length(which(all.pm.cov[[i]]>4))))
        allcells.names <- c(allcells.names,rep(allcells[i],length(which(all.pm.cov[[i]]>4))))
	allcells.types <- c(allcells.types,rep(all.celltypes[i],length(which(all.pm.cov[[i]]>4))))
	allcells.meth <- c(allcells.meth,all.pu.meth[[i]][which(all.pu.cov[[i]]>4)])
        allcells.regions <- c(allcells.regions,rep('PU',length(which(all.pu.cov[[i]]>4))))
        allcells.names <- c(allcells.names,rep(allcells[i],length(which(all.pu.cov[[i]]>4))))
	allcells.types <- c(allcells.types,rep(all.celltypes[i],length(which(all.pu.cov[[i]]>4))))
}

allcells.df <- data.frame(cell=allcells.names,region=allcells.regions,meth=allcells.meth,type=allcells.types)

library(ggplot2)

source('scMandTglobal.R')
classlookup <- cbind(colnames(full.countmatrix),dmr.clusters)
rownames(classlookup) <- classlookup[,1]
allcells.df$class <- classlookup[allcells.df$cell,2]

dm.means <- sapply(unique(allcells.df$cell),function(x)mean(allcells.df[which(allcells.df$cell==x & allcells.df$region=='DM'),'meth']))
cells.dmordered <- as.character(unique(allcells.df$cell)[order(dm.means)])

# need to downsample: take representative ranks from each set of cells (2I, serum)?
downsample.ranks <- round(quantile(1:length(dm.means),probs=c(0:11)/11))
downsample.cells.ordered <- as.character(unique(allcells.df$cell))[sapply(downsample.ranks,function(x)which(rank(dm.means)==x))]

# better. now try re-ordering factor(cell) based on mean

downsample.df <- allcells.df[which(as.character(allcells.df$cell) %in% downsample.cells.ordered),]
downsample.df$cell <- factor(downsample.df$cell,levels=downsample.cells.ordered)
downsample.df$Cluster <- apply(downsample.df,1,function(x)ifelse(x[4]=='2I','2i',ifelse(x[5]==1,'Naive-like','Primed-like')))

# now try 4 of each cluster
allcells.df$Cluster <- apply(allcells.df,1,function(x)ifelse(x[4]=='2I','2i',ifelse(x[5]==1,'Naive-like','Primed-like')))
downsample2.cells <- unlist(lapply(unique(allcells.df$Cluster),function(x)sample(unique(as.character(allcells.df$cell[which(allcells.df$Cluster==x)])),5)))
downsample2.means <- dm.means[sapply(downsample2.cells,function(x)which(as.character(unique(allcells.df$cell))==x))]
downsample2.cells.ordered <- downsample2.cells[order(downsample2.means)]
downsample2.df <- allcells.df[which(as.character(allcells.df$cell) %in% downsample2.cells.ordered),]
downsample2.df$cell <- factor(downsample2.df$cell,levels=downsample2.cells.ordered)

png(file='scMandT_enhanceletmeth_violinplot_DM_downsampled.png',width=12,height=6,units='in',res=300)
print(ggplot(downsample2.df[which(downsample2.df$region=='DM'),],aes(x=cell,y=meth)) + geom_violin(aes(fill=factor(Cluster)),scale='width') + scale_fill_manual(values=c('cyan','blue','orange')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()

png(file='scMandT_enhanceletmeth_violinplot_PU_downsampled.png',width=12,height=6,units='in',res=300)
print(ggplot(downsample2.df[which(downsample2.df$region=='PU'),],aes(x=cell,y=meth)) + geom_violin(aes(fill=factor(Cluster)),scale='width') + scale_fill_manual(values=c('cyan','blue','orange')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()

png(file='scMandT_enhanceletmeth_violinplot_PM_downsampled.png',width=12,height=6,units='in',res=300)
print(ggplot(downsample2.df[which(downsample2.df$region=='PM'),],aes(x=cell,y=meth)) + geom_violin(aes(fill=factor(Cluster)),scale='width') + scale_fill_manual(values=c('cyan','blue','orange')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()

