# load PM/DU/PU regions

source('scMandTglobal.R')

setwd('data/DNMT')

# obtain coverage files from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE75974&format=file 
# for speed of processing, CpG-level read count files were subsetted to regions of interest using 'bedtools intersect', resulting in a directory of files corresponding to each region class (DMR, PU or PM)
dnmt.dm.files <- list.files("DMRcov_mm10")
dnmt.dm.meth <- list()
dnmt.dm.cov <- list()
for(i in 1:length(dnmt.dm.files)){
        cat(paste("getting DM methylation for cell",i,"of",length(dnmt.dm.files),"\n"))
        this.bed <- read.table(paste("DMRcov_mm10/",dnmt.dm.files[i],sep=""))
        this.gr <- GRanges(seqnames=as.character(this.bed[[1]]),ranges=IRanges(start=this.bed[[2]],end=this.bed[[3]]),M=this.bed[[5]],U=this.bed[[6]])
        dnmt.dm.meth[[i]] <- sapply(DM.gr,getMethRegion,y=this.gr)
        dnmt.dm.cov[[i]] <- sapply(DM.gr,getCovRegion,y=this.gr)
}

dnmt.pu.files <- list.files("PUcov_mm10")
dnmt.pu.meth <- list()
dnmt.pu.cov <- list()
for(i in 1:length(dnmt.pu.files)){
        cat(paste("getting PU methylation for cell",i,"of",length(dnmt.pu.files),"\n"))
        this.bed <- read.table(paste("PUcov_mm10/",dnmt.pu.files[i],sep=""))
        this.gr <- GRanges(seqnames=as.character(this.bed[[1]]),ranges=IRanges(start=this.bed[[2]],end=this.bed[[3]]),M=this.bed[[5]],U=this.bed[[6]])
        dnmt.pu.meth[[i]] <- sapply(PU.gr,getMethRegion,y=this.gr)
        dnmt.pu.cov[[i]] <- sapply(PU.gr,getCovRegion,y=this.gr)
}

dnmt.pm.files <- list.files("PMcov_mm10")
dnmt.pm.meth <- list()
dnmt.pm.cov <- list()
for(i in 1:length(dnmt.pm.files)){
        cat(paste("getting PM methylation for cell",i,"of",length(dnmt.pm.files),"\n"))
        this.bed <- read.table(paste("PMcov_mm10/",dnmt.pm.files[i],sep=""))
        this.gr <- GRanges(seqnames=as.character(this.bed[[1]]),ranges=IRanges(start=this.bed[[2]],end=this.bed[[3]]),M=this.bed[[5]],U=this.bed[[6]])
        dnmt.pm.meth[[i]] <- sapply(PM.gr,getMethRegion,y=this.gr)
        dnmt.pm.cov[[i]] <- sapply(PM.gr,getCovRegion,y=this.gr)
}

# compute variance of methylation within DMs, PMs, PUs, for each cell
# filter out CpGs with too-low coverage (mark NA)

# compare 2i, serum, DNMT-ko, Tet-ko

alldm.cellCode = rep("serum",length(all.dm.files))
alldm.cellCode[grep("2I",all.dm.files)] <- "2i"
newdm.cellCode = rep('DNMT',length(dnmt.dm.files))
newdm.cellCode[grep('TET',dnmt.dm.files)] <- 'TET'
alldm.cellCode <- c(alldm.cellCode,newdm.cellCode)

old.dm.meth <- list()
old.dm.cov <- list()
for(i in 1:length(all.dm.gr)){
	cat(paste('computing DM methylation for cell',i,'of',length(all.dm.gr),'\n'))
	old.dm.meth[[i]] <- sapply(DM.gr,getMethRegion,y=all.dm.gr[[i]])
	old.dm.cov[[i]] <- sapply(DM.gr,getCovRegion,y=all.dm.gr[[i]])
}

old.pu.meth <- list()
old.pu.cov <- list()
for(i in 1:length(all.pu.gr)){
        cat(paste('computing PU methylation for cell',i,'of',length(all.pu.gr),'\n'))
        old.pu.meth[[i]] <- sapply(PU.gr,getMethRegion,y=all.pu.gr[[i]])
        old.pu.cov[[i]] <- sapply(PU.gr,getCovRegion,y=all.pu.gr[[i]])
}

old.pm.meth <- list()
old.pm.cov <- list()
for(i in 1:length(all.pm.gr)){
        cat(paste('computing PM methylation for cell',i,'of',length(all.pm.gr),'\n'))
        old.pm.meth[[i]] <- sapply(PM.gr,getMethRegion,y=all.pm.gr[[i]])
        old.pm.cov[[i]] <- sapply(PM.gr,getCovRegion,y=all.pm.gr[[i]])
}

olddm.var.min4 <- rep(NA,length(all.dm.files))
for(i in 1:length(all.dm.files)){
	olddm.var.min4[i] <- sd(old.dm.meth[[i]][old.dm.cov[[i]]>3])^2
}
newdm.var.min4 <- rep(NA,length(dnmt.dm.files))
for(i in 1:length(dnmt.dm.files)){
	newdm.var.min4[i] <- sd(dnmt.dm.meth[[i]][dnmt.dm.cov[[i]]>3])^2
}
oldpu.var.min4 <- rep(NA,length(all.pu.files))
for(i in 1:length(all.pu.files)){
        oldpu.var.min4[i] <- sd(old.pu.meth[[i]][old.pu.cov[[i]]>3])^2
}
newpu.var.min4 <- rep(NA,length(dnmt.pu.files))
for(i in 1:length(dnmt.pu.files)){
        newpu.var.min4[i] <- sd(dnmt.pu.meth[[i]][dnmt.pu.cov[[i]]>3])^2
}
oldpm.var.min4 <- rep(NA,length(all.pm.files))
for(i in 1:length(all.pm.files)){
        oldpm.var.min4[i] <- sd(old.pm.meth[[i]][old.pm.cov[[i]]>3])^2
}
newpm.var.min4 <- rep(NA,length(dnmt.pm.files))
for(i in 1:length(dnmt.pm.files)){
        newpm.var.min4[i] <- sd(dnmt.pm.meth[[i]][dnmt.pm.cov[[i]]>3])^2
}

allvar.df <- data.frame(variance=c(olddm.var.min4,newdm.var.min4,oldpu.var.min4,newpu.var.min4,oldpm.var.min4,newpm.var.min4),region=c(rep('DM',length(olddm.var.min4)+length(newdm.var.min4)),rep('PU',length(oldpu.var.min4)+length(newpu.var.min4)),rep('PM',length(oldpm.var.min4)+length(newpm.var.min4))),celltype=rep(alldm.cellCode,3))

library(ggplot2)

# repeat but with no background, and order PU/DM and 2i/serum/dnmtKO/tetKO
png(file='DNMT_VarianceDistributions_figure.png',width=12,height=6,units='in',res=300)
print(ggplot(allvar.df[which(allvar.df$region %in% c('PU','DM')),],aes(factor(celltype,level=c('2I','serum','DNMT','TET')),variance)) + geom_boxplot(aes(fill=factor(region,level=c('PU','DM')))) + scale_fill_manual(values=c('lightgreen','magenta')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()

# try subdividing serum cells
allvar.df$celltype2 <- as.character(allvar.df$celltype)
allvar.df$celltype2[which(allvar.df$celltype %in% c('2i','serum'))] <- rep(c('naive','primed')[dmr.clusters],3)
allvar.df$celltype2[which(allvar.df$celltype=='2i')] <- '2i'

# repeat but with no background, exclude knockouts: PU/DM/INT 2i/naive-like/primed-like

png(file='DNMT_VarianceDistributions_figure2.png',width=12,height=6,units='in',res=300)
print(ggplot(allvar.df[which(allvar.df$celltype2 %in% c('2I','naive','primed')),],aes(factor(celltype2,level=c('2I','naive','primed')),variance)) + geom_boxplot(aes(fill=factor(region,level=c('PU','DM','INT')))) + scale_fill_manual(values=c('lightgreen','magenta')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()

# now try for mean methylation, not variance

olddm.mean.min4 <- rep(NA,length(all.dm.files))
for(i in 1:length(all.dm.files)){
	        olddm.mean.min4[i] <- mean(old.dm.meth[[i]][old.dm.cov[[i]]>3])
}
newdm.mean.min4 <- rep(NA,length(dnmt.dm.files))
for(i in 1:length(dnmt.dm.files)){
	        newdm.mean.min4[i] <- mean(dnmt.dm.meth[[i]][dnmt.dm.cov[[i]]>3])
}
oldpu.mean.min4 <- rep(NA,length(all.pu.files))
for(i in 1:length(all.pu.files)){
	        oldpu.mean.min4[i] <- mean(old.pu.meth[[i]][old.pu.cov[[i]]>3])
}
newpu.mean.min4 <- rep(NA,length(dnmt.pu.files))
for(i in 1:length(dnmt.pu.files)){
	        newpu.mean.min4[i] <- mean(dnmt.pu.meth[[i]][dnmt.pu.cov[[i]]>3])
}
oldpm.mean.min4 <- rep(NA,length(all.pm.files))
for(i in 1:length(all.pm.files)){
	        oldpm.mean.min4[i] <- mean(old.pm.meth[[i]][old.pm.cov[[i]]>3])
}
newpm.mean.min4 <- rep(NA,length(dnmt.pm.files))
for(i in 1:length(dnmt.pm.files)){
	        newpm.mean.min4[i] <- mean(dnmt.pm.meth[[i]][dnmt.pm.cov[[i]]>3])
}

allmean.df <- data.frame(meanmeth=c(olddm.mean.min4,newdm.mean.min4,oldpu.mean.min4,newpu.mean.min4,oldpm.mean.min4,newpm.mean.min4),region=c(rep('DM',length(olddm.mean.min4)+length(newdm.mean.min4)),rep('PU',length(oldpu.mean.min4)+length(newpu.mean.min4)),rep('PM',length(oldpm.mean.min4)+length(newpm.mean.min4))),celltype=rep(alldm.cellCode,3))
allmean.df$celltype2 <- as.character(allmean.df$celltype)
allmean.df$celltype2[which(allmean.df$celltype %in% c('2i','serum'))] <- rep(c('naive','primed')[dmr.clusters],3)
allmean.df$celltype2[which(allmean.df$celltype=='2i')] <- '2i'

png(file='DNMT_MeanMethDistributions_figure.png',width=12,height=6,units='in',res=300)
print(ggplot(allmean.df[which(allmean.df$region %in% c('PU','DM')),],aes(factor(celltype,level=c('2I','serum','DNMT','TET')),variance)) + geom_boxplot(aes(fill=factor(region,level=c('PU','DM')))) + scale_fill_manual(values=c('lightgreen','magenta')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()

png(file='DNMT_MeanMethDistributions_figure2.png',width=12,height=6,units='in',res=300)
print(ggplot(allmean.df[which(allmean.df$celltype2 %in% c('2I','naive','primed')),],aes(factor(celltype2,level=c('2I','naive','primed')),variance)) + geom_boxplot(aes(fill=factor(region,level=c('PU','DM','INT')))) + scale_fill_manual(values=c('lightgreen','magenta')) + theme(panel.background=element_rect(fill='white',colour='black')))
dev.off()

