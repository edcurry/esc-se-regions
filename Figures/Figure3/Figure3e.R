# read ChIPseq files table

library(ShortRead)

chip.files <- read.table("allChIPseqFiles.txt",sep=",",head=T)
chip.files <- rbind(chip.files,read.table("allChIPseqFiles2.txt",sep=",",head=T))
chip.files <- chip.files[!is.na(as.character(chip.files$File)),]
chip.files <- chip.files[which(chip.files$ChIP %in% c('ATAC','Pou5f1','H3K27ac')),]
chip.files <- chip.files[!duplicated(paste0(chip.files$Cell,chip.files$ChIP)),]

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
	this.ChIP.logFC[which(this.ChIP.rpkm<((1000*4*1e6/length(this.ChIP.reads))/width(enhancelet.gr)) & this.input.rpkm<((1000*4*1e6/length(this.input.reads))/width(enhancelet.gr)))] <- NA
	enhancelet.matrix[,i] <- this.ChIP.logFC
	}
}

# load methylation data from suppl table S3
# make dataframe with Cell (ESC:'pct_meth_ser-ESC_Marks', EpiLC-activin:'pct_meth_EpiLC_Surani', EpiSC:'pct_meth_EpiSC_Veillard'), 'type' (PU, DM)

# make dataframe with Cell (ESC, EpiLC-activin, EpiSC), featureType (PU, DM) 
ATACdf <- data.frame(ATAC.logFC=as.numeric(enhancelet.matrix[,which(chip.files$ChIP=='ATAC')]),Cell=rep(chip.files[which(chip.files$ChIP=='ATAC'),'Cell'],each=nrow(enhancelet.matrix)),region=featureType)
library(ggplot2)
png(file='Fig3e_ATAC.png',height=8,width=8,units='in',res=300)
print(ggplot(ATACdf[which(ATACdf$region %in% c('PU','DM')),],aes(x=factor(Cell,level=c('ESC','EpiLC','EpiSC')),y=ATAC.logFC)) + geom_boxplot(aes(fill=factor(region,level=c('PU','DM')))) + scale_fill_manual(values=c('lightgreen','magenta')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()

Pou5f1df <- data.frame(Pou5f1.logFC=as.numeric(enhancelet.matrix[,which(chip.files$ChIP=='Pou5f1')]),Cell=rep(chip.files[which(chip.files$ChIP=='Pou5f1'),'Cell'],each=nrow(enhancelet.matrix)),region=featureType)
png(file='Fig3e_Pou5f1.png',height=8,width=8,units='in',res=300)
print(ggplot(Pou5f1df[which(Pou5f1df$region %in% c('PU','DM')),],aes(x=factor(Cell,level=c('ESC','EpiLC','EpiSC')),y=Pou5f1.logFC)) + geom_boxplot(aes(fill=factor(region,level=c('PU','DM')))) + scale_fill_manual(values=c('lightgreen','magenta')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()

H3K27acdf <- data.frame(H3K27ac.logFC=as.numeric(enhancelet.matrix[,which(chip.files$ChIP=='H3K27ac')]),Cell=rep(chip.files[which(chip.files$ChIP=='H3K27ac'),'Cell'],each=nrow(enhancelet.matrix)),region=featureType)
png(file='Fig3e_H3K27ac.png',height=8,width=8,units='in',res=300)
print(ggplot(H3K27acdf[which(H3K27acdf$region %in% c('PU','DM')),],aes(x=factor(Cell,level=c('ESC','EpiLC','EpiSC')),y=H3K27ac.logFC)) + geom_boxplot(aes(fill=factor(region,level=c('PU','DM')))) + scale_fill_manual(values=c('lightgreen','magenta')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()


