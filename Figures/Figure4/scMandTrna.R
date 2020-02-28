source("scMandTglobal.R")

# find differentially expressed genes according to dmr cluster grouping

# RNAseq counts table from: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74534&format=file&file=GSE74534%5FRNA%2Dseq%5Fnormalized%5Fcounts%2Etxt%2Egz
all.RNAseq <- as.matrix(read.table("GSE74534_RNA-seq_normalized_counts.txt",sep="\t",head=T,row.names=1))

RNAseq.expressed <- all.RNAseq[which(rowSums(all.RNAseq>1)>5),]
RNAseq.expressed <- log(RNAseq.expressed,base=2)
RNAseq.expressed[RNAseq.expressed==(-Inf)] <- min(RNAseq.expressed[RNAseq.expressed>(-Inf)],na.rm=T)
colnames(RNAseq.expressed) = sapply(colnames(RNAseq.expressed),function(x)paste(strsplit(x,split="_")[[1]][2:3],collapse="_"))

candidategenes <- c("Esrrb","Nr5a2","Klf2","Klf4","Pou5f1","Nanog")

library(biomaRt)
mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",host="feb2014.archive.ensembl.org")
refgenes <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end'), mart = mouse)
refgenes <- refgenes[!duplicated(refgenes$ensembl_gene_id),]
rownames(refgenes) <- as.character(refgenes$ensembl_gene_id)

RNAseq.symbols <- refgenes[rownames(RNAseq.expressed),"mgi_symbol"]

# make boxplots of candidates (separate naive-like and primed-like)

# plot distribution of marker genes in 2i, naive-like and primed-like subsets
library(ggplot2)
names(dmr.clusters) <- toupper(names(dmr.clusters))
candidatemarker.ensIDs <- rownames(RNAseq.expressed)[which(RNAseq.symbols %in% candidategenes)]
for(i in 1:length(candidatemarker.ensIDs)){
        boxplotdf <- data.frame(cell=colnames(RNAseq.expressed),subtype=c('naive-like','primed-like')[dmr.clusters[colnames(RNAseq.expressed)]],expression=RNAseq.expressed[candidatemarker.ensIDs[i],])
#       colnames(boxplotdf)[3] <- paste(RNAseq.symbols[which(rownames(RNAseq.expressed)==candidatemarker.ensIDs[i])],'_expression',sep='')
        png(file=paste('scMandT_ExpressionBySubclass_',candidatemarker.ensIDs[i],'.png',sep=''),width=5,height=5,units='in',res=300)
# add specific colours for fill aesthetic (naive-like=blue, primed-like=orange)
        print(ggplot(boxplotdf,aes(factor(subtype),expression)) + geom_boxplot(aes(fill=factor(subtype))) + scale_fill_manual(values=c('blue','orange')) + ggtitle(RNAseq.symbols[which(rownames(RNAseq.expressed)==candidatemarker.ensIDs[i])]) + theme(panel.background=element_rect(fill='white',colour='black')))
        dev.off()
}

# assess differential expression

library(limma)
design <- cbind(intercept=1,dm.class=dmr.clusters-1)
rownames(design) <- toupper(rownames(design))
design <- design[colnames(RNAseq.expressed),]
##v <- voom(RNAseq.expressed,design=design)
vfit <- lmFit(RNAseq.expressed,design)
efit <- eBayes(vfit)
dmrclass.diffexp <- topTable(efit,coef=2,number=nrow(RNAseq.expressed))
dmrclass.diffexp$symbol = refgenes[rownames(dmrclass.diffexp),"mgi_symbol"]

genes.RNAseq <- refgenes[rownames(RNAseq.expressed),"mgi_symbol"]

naivepluri <- as.character(read.table("~/Emma/Smiths/naivePluriGenes.txt")[[1]])
naivepluri <- intersect(naivepluri,as.character(genes.RNAseq))
generalpluri <- as.character(read.table("~/Emma/Smiths/generalPluriGenes.txt")[[1]])
generalpluri <- intersect(generalpluri,as.character(genes.RNAseq))
postimplant <- as.character(read.table("~/Emma/Smiths/postImplantGenes.txt")[[1]])
postimplant <- intersect(postimplant,as.character(genes.RNAseq))

plotmat3 <- t(RNAseq.expressed[sapply(c(naivepluri,generalpluri,postimplant),function(x)which(genes.RNAseq==x)),])
colnames(plotmat3) <- c(naivepluri,generalpluri,postimplant)

png(file="scRNAseqGXsignatureROC.png",width=6,height=6,units='in',res=300)
plot(x=c(0,cumsum(design[order(colMeans(p3.z[naivepluri,]),decreasing=T),2]))/sum(design[,2]),y=c(0,cumsum(1-design[order(colMeans(p3.z[naivepluri,]),decreasing=T),2]))/sum(1-design[,2]),type="l",lwd=2,col="yellow",xlim=c(0,1),ylim=c(0,1),ylab="sensitivity",xlab="1-specificity")
points(x=c(0,cumsum(design[order(colMeans(p3.z[generalpluri,]),decreasing=T),2]))/sum(design[,2]),y=c(0,cumsum(1-design[order(colMeans(p3.z[generalpluri,]),decreasing=T),2]))/sum(1-design[,2]),type="l",lwd=2,col="orange")
points(x=c(0,cumsum(design[order(colMeans(p3.z[postimplant,]),decreasing=T),2]))/sum(design[,2]),y=c(0,cumsum(1-design[order(colMeans(p3.z[postimplant,]),decreasing=T),2]))/sum(1-design[,2]),type="l",lwd=2,col="red")
legend("bottomright",legend=c("naive-pluri","general-pluri","post-implant"),lwd=2,col=c("yellow","orange","red"))
dev.off()

serum <- grep("SERUM",names(dmr.clusters))

png(file="scRNAseqGXcandidatesROC2.png",width=6,height=6,units='in',res=300)
plot(x=c(0,cumsum((dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Nr5a2"),],decreasing=T)]))/sum((dmr.clusters[colnames(RNAseq.expressed)]-1)),y=c(0,cumsum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Nr5a2"),],decreasing=T)]))/sum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)),type="l",lwd=2,col="cyan",xlim=c(0,1),ylim=c(0,1),ylab="sensitivity",xlab="1-specificity")
points(x=c(0,cumsum((dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Esrrb"),],decreasing=T)]))/sum((dmr.clusters[colnames(RNAseq.expressed)]-1)),y=c(0,cumsum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Esrrb"),],decreasing=T)]))/sum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)),type="l",lwd=2,col="blue")
points(x=c(0,cumsum((dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Klf2"),],decreasing=T)]))/sum((dmr.clusters[colnames(RNAseq.expressed)]-1)),y=c(0,cumsum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Klf2"),],decreasing=T)]))/sum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)),type="l",lwd=2,col="yellow")
points(x=c(0,cumsum((dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Klf4"),],decreasing=T)]))/sum((dmr.clusters[colnames(RNAseq.expressed)]-1)),y=c(0,cumsum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Klf4"),],decreasing=T)]))/sum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)),type="l",lwd=2,col="green")
points(x=c(0,cumsum((dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Nanog"),],decreasing=T)]))/sum((dmr.clusters[colnames(RNAseq.expressed)]-1)),y=c(0,cumsum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Nanog"),],decreasing=T)]))/sum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)),type="l",lwd=2,col="red")
points(x=c(0,cumsum((dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Pou5f1"),],decreasing=T)]))/sum((dmr.clusters[colnames(RNAseq.expressed)]-1)),y=c(0,cumsum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Pou5f1"),],decreasing=T)]))/sum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)),type="l",lwd=2,col="purple")
points(x=c(0,cumsum((dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Zfp42"),],decreasing=T)]))/sum((dmr.clusters[colnames(RNAseq.expressed)]-1)),y=c(0,cumsum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)[order(RNAseq.expressed[which(RNAseq.symbols=="Zfp42"),],decreasing=T)]))/sum(1-(dmr.clusters[colnames(RNAseq.expressed)]-1)),type="l",lwd=2,col="magenta")
legend("bottomright",legend=c("Esrrb","Klf2","Klf4","Nanog","Nr5a2","Pou5f1","Zfp42"),lwd=2,col=c("blue","yellow","green","red","cyan","purple","magenta"))
dev.off()

all.fcs <- apply(plotmat3,2,function(x)median(x[design[,2]==0])-median(x[design[,2]==1]))

makeROC <- function(gene){
	        out <- NULL
        if(length(gene)>1){
		                out <- data.frame(y=c(0,cumsum(1-design[order(colMeans(p3.z[gene,]),decreasing=T),2]))/sum(1-design[,2]),x=c(0,cumsum(design[order(colMeans(p3.z[gene,]),decreasing=T),2])/sum(design[,2])))
	        }
	        if(length(gene)==1){
			                out <- data.frame(y=c(0,cumsum(1-design[order(p3.z[gene,],decreasing=T),2]))/sum(1-design[,2]),x=c(0,cumsum(design[order(p3.z[gene,],decreasing=T),2])/sum(design[,2])))
		        }
		        out
}

getAUC <- function(ROC){
	        sum((1-ROC$x[1:(nrow(ROC)-1)])*(ROC$y[2:nrow(ROC)]-ROC$y[1:(nrow(ROC)-1)]))
}

# make table S4b
tabs4b.roc <- lapply(list(Naive=naivepluri,General=generalpluri,Primed=postimplant,"Esrrb","Klf2","Klf4","Nanog","Pou5f1","Zfp42"),makeROC)
tabs4b.auc <- sapply(tabs4b.roc,getAUC)

# include permutation test for p-values
randomAUCs <- rep(NA,1e5)
for(i in 1:length(randomAUCs)){
	        randomOrder <- sample(ncol(p3.z))
        randomAUCs[i] <- getAUC(data.frame(y=c(0,cumsum(1-design[randomOrder,2]))/sum(1-design[,2]),x=c(0,cumsum(design[randomOrder,2])/sum(design[,2]))))
}
tabs4b.pvals <- (sapply(tabs4b.auc,function(x)sum(randomAUCs>=x))+1)/(1e5+1)

tabs4b <- data.frame(Gene.set=c("Naive","General","Primed","Esrrb","Klf2","Klf4","Nanog","Pou5f1","Zfp42"),AUC=tabs4b.auc,p.value=tabs4b.pvals)
write.table(tabs4b,file="TableS4b_GeneAUCs.txt",quote=F,sep="\t",row.names=F)

p4.z <- (RNAseq.expressed[which(genes.RNAseq!=""),]-apply(RNAseq.expressed[which(genes.RNAseq!=""),],1,median,na.rm=T))/apply(RNAseq.expressed[which(genes.RNAseq!=""),],1,sd,na.rm=T)

fullauc <- function(x){
	        ROC <- data.frame(y=c(0,cumsum(1-design[order(x,decreasing=T),2]))/sum(1-design[,2]),x=c(0,cumsum(design[order(x,decreasing=T),2])/sum(design[,2])))
        sum((1-ROC$x[1:(nrow(ROC)-1)])*(ROC$y[2:nrow(ROC)]-ROC$y[1:(nrow(ROC)-1)]))
}
genomewide.aucs <- apply(p4.z,1,fullauc)
genomewide.auc.df <- data.frame(Gene=genes.RNAseq[which(genes.RNAseq!="")],AUC=genomewide.aucs,FC=2^all.fcs[which(genes.RNAseq!="")])
genomewide.auc.df <- genomewide.auc.df[order(genomewide.auc.df$AUC,decreasing=T),]
write.table(genomewide.auc.df,file="scMandT_2ilikeROCAUC_genomewide.txt",sep="\t",row.names=F,quote=F)

class1genes <- c("Klf13","Lefty1","Med13l","Pou5f1","Nanog","Smarcad1","Tet1","Sox2","Otx2")
class2genes <- c("Esrrb","Klf2","Klf4","Klf5","Tet2","Tdh","Tbx3","Tfcp2l1","Zfp42")

dmrclass.gx.clust <- hclust(dist(t(RNAseq.expressed[rownames(dmrclass.diffexp)[1:20],])))
library(NMF)
rnaseq.plotmat <- t(RNAseq.expressed[c(which(genes.RNAseq %in% class1genes),which(genes.RNAseq %in% class2genes)),])
colnames(rnaseq.plotmat) <- c(genes.RNAseq[which(genes.RNAseq %in% class1genes)],genes.RNAseq[which(genes.RNAseq %in% class2genes)])
rnaseq.plotmat <- rnaseq.plotmat[,c(sort(colnames(rnaseq.plotmat)[colnames(rnaseq.plotmat) %in% class1genes]),sort(colnames(rnaseq.plotmat)[colnames(rnaseq.plotmat) %in% class2genes]))]
png(file="scMandT_2iLike_ClassIvIIexpressionHeatmap2.png",width=10,height=8,units='in',res=200)
aheatmap(rnaseq.plotmat,col="-RdBu:100",annRow=list(state=c("2i-like","other")[1+design[,2]]),labRow=NA,Colv=NA,Rowv=as.dendrogram(dmrclass.gx.clust),scale="col",revC=T,annColors=list(state=c("blue","orange")))
dev.off()

png(file="scMandT_2iLike_ClassIvIIexpressionBoxplot2.png",width=8,height=8,units='in',res=200)
par(mfrow=c(1,2))
boxplot(list(class1=(RNAseq.expressed[which(genes.RNAseq %in% class1genes),design[,2]==0]-apply(RNAseq.expressed[which(genes.RNAseq %in% class1genes),],MARGIN=1,median))/apply(RNAseq.expressed[which(genes.RNAseq %in% class1genes),],MARGIN=1,sd),class2=(RNAseq.expressed[which(genes.RNAseq %in% class2genes),design[,2]==0]-apply(RNAseq.expressed[which(genes.RNAseq %in% class2genes),],MARGIN=1,median))/apply(RNAseq.expressed[which(genes.RNAseq %in% class2genes),],MARGIN=1,sd)),ylab="gene expression z-score",col=c("lightgreen","pink"),main="2i-like")
boxplot(list(class1=(RNAseq.expressed[which(genes.RNAseq %in% class1genes),design[,2]==1]-apply(RNAseq.expressed[which(genes.RNAseq %in% class1genes),],MARGIN=1,median))/apply(RNAseq.expressed[which(genes.RNAseq %in% class1genes),],MARGIN=1,sd),class2=(RNAseq.expressed[which(genes.RNAseq %in% class2genes),design[,2]==1]-apply(RNAseq.expressed[which(genes.RNAseq %in% class2genes),],MARGIN=1,median))/apply(RNAseq.expressed[which(genes.RNAseq %in% class2genes),],MARGIN=1,sd)),ylab="gene expression z-score",col=c("lightgreen","pink"),main="other")
dev.off()


