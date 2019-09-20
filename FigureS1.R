library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm9, quietly = TRUE)

esc.se <- read.table(file="data/ESC_SE_locations_trimmed.txt",sep="\t",comment.char="",quote="",head=F,stringsAsFactors=F)
colnames(esc.se) <- c("chr","start","stop")
se.ranges <- with(esc.se, GRanges(chr, IRanges(start=start, end=stop)))

se.seqs <- getSeq(Mmusculus, seqnames(se.ranges), start(se.ranges), end(se.ranges), as.character = T)
se.length <- width(se.ranges)

esc.te <- read.table(file="data/TE_locations_trimmed.txt",sep="\t",comment.char="",quote="",head=F,stringsAsFactors=F)
colnames(esc.te) <- c("chr","start","stop")
te.ranges <- with(esc.te, GRanges(chr, IRanges(start=start, end=stop)))

te.seqs <- getSeq(Mmusculus, seqnames(te.ranges), start(te.ranges), end(te.ranges), as.character = T)
te.length <- width(te.ranges)

prob.se <- read.table(file="data/ProB_SE_locations_trimmed.txt",sep="\t",comment.char="",quote="",head=F,stringsAsFactors=F)
colnames(prob.se) <- c("chr","start","stop")
prob.ranges <- with(prob.se, GRanges(chr, IRanges(start=start, end=stop)))

prob.seqs <- getSeq(Mmusculus, seqnames(prob.ranges), start(prob.ranges), end(prob.ranges), as.character = T)
prob.length <- width(prob.ranges)

#calculate CpG, C and G frequency
library(stringr)
se.cg.count <- rep(NA,nrow(esc.se))
se.c.count <- rep(NA,nrow(esc.se))
se.g.count <- rep(NA,nrow(esc.se))

for(i in 1:nrow(esc.se)){
se.cg.count[i] <- str_count(se.seqs[i], fixed("CG"))+str_count(se.seqs[i], fixed("GC"))
se.c.count[i] <- str_count(se.seqs[i], fixed("C"))
se.g.count[i] <- str_count(se.seqs[i], fixed("G"))
}

se.CpG <- se.cg.count/se.length
se.GC <- (se.c.count+se.g.count)/se.length

te.cg.count <- rep(NA,nrow(esc.te))
te.c.count <- rep(NA,nrow(esc.te))
te.g.count <- rep(NA,nrow(esc.te))

for(i in 1:nrow(esc.te)){
te.cg.count[i] <- str_count(te.seqs[i], fixed("CG"))+str_count(te.seqs[i], fixed("GC"))
te.c.count[i] <- str_count(te.seqs[i], fixed("C"))
te.g.count[i] <- str_count(te.seqs[i], fixed("G"))
}

te.CpG <- te.cg.count/te.length
te.GC <- (te.c.count+te.g.count)/te.length

prob.cg.count <- rep(NA,nrow(prob.se))
prob.c.count <- rep(NA,nrow(prob.se))
prob.g.count <- rep(NA,nrow(prob.se))

for(i in 1:nrow(prob.se)){
prob.cg.count[i] <- str_count(prob.seqs[i], fixed("CG"))+str_count(prob.seqs[i], fixed("GC"))
prob.c.count[i] <- str_count(prob.seqs[i], fixed("C"))
prob.g.count[i] <- str_count(prob.seqs[i], fixed("G"))
}

prob.CpG <- prob.cg.count/prob.length
prob.GC <- (prob.c.count+prob.g.count)/prob.length

# make plots

png(file="SupplFigS1_length.png",width=8,height=8,units='in',res=300)
boxplot(list(ESC.SE=se.length/1000,ESC.TE=te.length/1000,ProB.SE=prob.length/1000),col=c('red','yellow','blue'),notch=T,ylab='Length (kb)',xlab='')
dev.off()

png(file="SupplFigS1_GC.png",width=8,height=8,units='in',res=300)
boxplot(list(ESC.SE=se.GC*100,ESC.TE=te.GC*100,ProB.SE=prob.GC*100),col=c('red','yellow','blue'),notch=T,ylab='%GC',xlab='')
dev.off()

png(file="SupplFigS1_CpG.png",width=8,height=8,units='in',res=300)
boxplot(list(ESC.SE=se.CpG,ESC.TE=te.CpG,ProB.SE=prob.CpG),col=c('red','yellow','blue'),notch=T,ylab='CpG density',xlab='',ylim=c(0,0.1))
dev.off()
