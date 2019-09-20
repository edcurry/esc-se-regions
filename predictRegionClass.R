# read ChIPseq files table

library(ShortRead)

# only use Dnmts, Tets and H3K4me3 ChIP
chip.files <- read.table("predictionChIPseqFiles.txt",sep=",",head=T)
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

# add CpG density to enhancelet.matrix

source('CpG_count.R')

enhancelet.df <- as.data.frame(enhancelet.matrix)
enhancelet.df$GCcontent <- enhancelets.GC
enhancelet.df$CpGdensity <- enhancelets.cg.count/(enhancelets.length) 
enhancelet.df$DM <- as.numeric(featureType=='DM')

# predict 'featureType' based on columns of enhancelet.matrix

basiclm.ESC <- glm(DM~.,data=enhancelet.df[featureType %in% c('DM','PU'),c(1:7,12:14)],family='binomial')
summary(basiclm.ESC)
#Coefficients:
#              Estimate Std. Error z value Pr(>|z|)
#(Intercept)     4.8089     1.5990   3.007 0.002634 **
#`ESC-Tet1`     -1.0605     0.3101  -3.420 0.000626 ***
#`ESC-Dnmt3a`    0.9074     0.6839   1.327 0.184609
#`ESC-Dnmt3b`    0.8289     0.7248   1.144 0.252787
#`ESC-H3K4me1`   1.0510     0.3521   2.985 0.002835 **
#`ESC-H3K4me3`  -0.3998     0.1433  -2.789 0.005284 **
#`ESC-H3K27ac`   0.5537     0.2275   2.434 0.014952 *
#`ESC-ATAC`     -0.9417     0.1431  -6.581 4.68e-11 ***
#GCcontent      -0.9116     3.5452  -0.257 0.797071
#CpGdensity    -33.8482    10.7502  -3.149 0.001641 **

mean(1-abs(enhancelet.df[featureType %in% c('DM','PU'),'DM']-as.numeric(predict(basiclm.ESC,type='response')>0.5)))
#[1] 0.8392857

# now try train/test split, using ESC data only (it appears DNMT is not important, nor is GCcontent)
enhancelet.ESdf <- enhancelet.df[featureType %in% c('DM','PU'),c(1:7,12:14)]
kfoldCV.acc <- rep(NA,10)
for(i in 1:10){
	this.test <- sample(nrow(enhancelet.ESdf),round(nrow(enhancelet.ESdf)/3))
	this.train <- setdiff(1:nrow(enhancelet.ESdf),this.test)
	this.esfit <- glm(DM~.,data=enhancelet.ESdf[this.train,],family='binomial')
	kfoldCV.acc[i] <- mean(1-abs(enhancelet.ESdf[this.test,'DM']-as.numeric(predict(this.esfit,newdata=enhancelet.ESdf[this.test,],type='response')>0.5)))
}
quantile(kfoldCV.acc,probs=c(0.05,0.5,0.95))
#       5%       50%       95% 
#0.8120536 0.8325893 0.8531250

# so, predicts with 83% accuracy (95%CI:0.81-0.85)
