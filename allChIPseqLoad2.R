# read ChIPseq files table

library(ShortRead)

chip.files <- read.table("allChIPseqFiles2.txt",sep=",",head=T)
chip.files <- chip.files[!is.na(as.character(chip.files$File)),]

# define regions

DM <- read.table("data/DMR_mm9.bed",sep="\t",head=F)
DM.gr <- GRanges(seqnames=as.character(DM[[1]]),ranges=IRanges(start=DM[[2]],end=DM[[3]]))
PM <- read.table("data/PM_mm9.bed",sep="\t",head=F)
PM.gr <- GRanges(seqnames=as.character(PM[[1]]),ranges=IRanges(start=PM[[2]],end=PM[[3]]))
PU <- read.table("data/PU_mm9.bed",sep="\t",head=F)
PU.gr <- GRanges(seqnames=as.character(PU[[1]]),ranges=IRanges(start=PU[[2]],end=PU[[3]]))
flankingUp <- read.table("data/flankingUp.bed",sep="\t",head=F)
flankingUp.gr <- GRanges(seqnames=as.character(flankingUp[[1]]),ranges=IRanges(start=floor(flankingUp[[2]]),end=ceiling(flankingUp[[3]])))
flankingDown <- read.table("data/flankingDown.bed",sep="\t",head=F)
flankingDown.gr <- GRanges(seqnames=as.character(flankingDown[[1]]),ranges=IRanges(start=floor(flankingDown[[2]]),end=ceiling(flankingDown[[3]])))
enhancelet.gr <- c(DM.gr,PM.gr,PU.gr,flankingUp.gr,flankingDown.gr)
featureType <- c(rep("DM",length(DM.gr)),rep("PM",length(PM.gr)),rep("PU",length(PU.gr)),rep("flU",length(flankingUp.gr)),rep("flD",length(flankingDown.gr)))

enhancelet.matrix <- array(NA,dim=c(length(enhancelet.gr),nrow(chip.files)))
colnames(enhancelet.matrix) <- paste(as.character(chip.files$Cell),"-",as.character(chip.files$ChIP),sep="")

# for each ChIP file, compute rpkm, normalize to input, add to matrix

for(i in 1:nrow(chip.files)){
	if(!is.na(as.character(chip.files$File)[i])){
	cat(paste("generating RPKMs for sample ",i,": ",as.character(chip.files$Cell)[i],"-",as.character(chip.files$ChIP)[i],"\n",sep=""))
	cat("reading ChIP bam \n")
	this.ChIP.reads <- readGAlignments(file=as.character(chip.files$File)[i],param=ScanBamParam())
	cat("counting overlaps \n")
	this.ChIP.counts <- summarizeOverlaps(features=enhancelet.gr,reads=this.ChIP.reads)#,mode="IntersectionNotEmpty")
	this.ChIP.rpkm <- (1000*assay(this.ChIP.counts)*1e6/length(this.ChIP.reads))/width(enhancelet.gr)
	this.ChIP.logFC <- log(this.ChIP.rpkm+1,base=2)
	
	if(!is.na(as.character(chip.files$Input)[i])){
		cat("reading input bam \n")
		this.input.reads <- readGAlignments(file=as.character(chip.files$Input)[i],param=ScanBamParam())
	        cat("counting overlaps \n")
		this.input.counts <- summarizeOverlaps(features=enhancelet.gr,reads=this.input.reads)#,mode="IntersectionNotEmpty")
		this.input.rpkm <- (1000*assay(this.input.counts)*1e6/length(this.input.reads))/width(enhancelet.gr)
		this.ChIP.logFC <- log(((this.ChIP.rpkm+1)/(this.input.rpkm+1)),base=2)
	}
	
	enhancelet.matrix[,i] <- this.ChIP.logFC
	}
}

# compute average enrichments per region-type

feature.matrix <- array(NA,dim=c(4,nrow(chip.files)))
rownames(feature.matrix) <- c("PU","DM","PM","Flank")
colnames(feature.matrix) <- paste(as.character(chip.files$Cell),"-",as.character(chip.files$ChIP),sep="")
feature.matrix[1,] <- apply(enhancelet.matrix[which(featureType=="PU"),],2,median,na.rm=T)
feature.matrix[2,] <- apply(enhancelet.matrix[which(featureType=="DM"),],2,median,na.rm=T)
feature.matrix[3,] <- apply(enhancelet.matrix[which(featureType=="PM"),],2,median,na.rm=T)
feature.matrix[4,] <- apply(enhancelet.matrix[which(featureType %in% c("flU","flD")),],2,median,na.rm=T)
feature.matrix.z <- t((t(feature.matrix)-apply(t(feature.matrix),1,median,na.rm=T))/apply(t(feature.matrix),1,sd,na.rm=T))

#write.table(feature.matrix,file="AllChIPseqFeatureMatrix_Epi.txt",sep="\t",col.names=NA,quote=F)

