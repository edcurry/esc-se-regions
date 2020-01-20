setwd("C:/Users/LucJouneau/Documents/INRA/Alice/Article Azuara/article/Figure1")

Klf4_SE_coordinates=c(4,55469259,55491081)
Klf13_SE_coordinates=c(7,71092246,71102481)
Lefty1_SE_coordinates=c(1,182854521,182864307)

INT_regions=read.table(file="data/INT_mm9.bed",header=FALSE,sep="\t")
DM_regions=read.table(file="data/DM_mm9.bed",header=FALSE,sep="\t")
PU_regions=read.table(file="data/PU_mm9.bed",header=FALSE,sep="\t")

#-----------------
# Methylation call
#-----------------

for (gene_se in c("Klf4","Klf13","Lefty1")) {

	maxGap=999

	if (gene_se=="Klf4") {SE_coordinates=Klf4_SE_coordinates}
	if (gene_se=="Klf13") {SE_coordinates=Klf13_SE_coordinates}
	if (gene_se=="Lefty1") {SE_coordinates=Lefty1_SE_coordinates}

	#-----------------
	# Methylation call
	#-----------------
	chr=SE_coordinates[1]
	for (cell_line in c("ESC","EpiSC")) {
		if (cell_line=="ESC") {
			tab=read.table(file=paste("data/cnv384_CpG_methcounts_chr",chr,".bed",sep=""),sep="\t",header=F)
		} else {
			tab=read.table(file=paste("data/29_EpiSC_8.73_merged.rmDup.methcounts_chr",chr,".bed",sep=""),sep="\t",header=F)
		}
	
		coverage=as.numeric(gsub("CpG:","",tab[,4]))
		keep=coverage>=7
		tabWithCoverage=tab[keep,]
		tabWithCoverage2tab=which(keep)
	
		distance=tabWithCoverage[2:nrow(tabWithCoverage),3]-tabWithCoverage[1:(nrow(tabWithCoverage)-1),2]
		start=c(1,(which(distance>maxGap)+1))
		end=start[2:length(start)]-1
		#On ajoute la derni?re coordonn?e de fin
		end=c(end,nrow(tabWithCoverage))
	
		#On conserve les r?gions qui ont au moins 10 positions
		segmentLength=end-start+1
		keepSegment=segmentLength>=10
		start=start[keepSegment]
		end=end[keepSegment]
		coords=rbind(start,end)
		obs=apply(coords,2,function(x) {tabWithCoverage[x[1]:x[2],5]})
	
		library(tileHMM)
		hmm.init <- hmm.setup(
			unlist(obs),
			state = c("non-methylated","methylated"),
			pos.state = 1
		)
		hmm.opt <- viterbiEM(hmm.init, obs)
		post <- lapply(obs, posterior, hmm.opt)
		state.seq <- lapply(post, apply, 2, which.max)
		state.seq <- states(hmm.opt)[c(state.seq, recursive=TRUE)]
	
		hmm_result=apply(coords,2,function(x) {tabWithCoverage[x[1]:x[2],]})
		hmm_result=do.call(rbind.data.frame,hmm_result)
		colnames(hmm_result)=colnames(tabWithCoverage)
		hmm_result=cbind(hmm_result,state.seq)
		colnames(hmm_result)=c("Chromosome","Start","End","Coverage","Methylation","Strand","State")
		hmm_result=hmm_result[,-6]
		if (cell_line=="ESC") {
			ESC_result=hmm_result
		} else {#Assert EpiSC cell line
			EpiSC_result=hmm_result
		}
	}

	#-----------------
	# Figure 2
	#-----------------
	jpeg(paste("Figure2a - ",gene_se,".jpg",sep=""))
	plot(x=c(),y=c(),
	     xlim=c(SE_coordinates[2],SE_coordinates[3]),ylim=c(1,11),
	     main=gene_se,
	     xlab="",ylab="",xaxt="n",yaxt="n",frame.plot=FALSE
	)
	

	for (cell_line in c("ESC","EpiSC")) {
		if( cell_line=="ESC") {
			hmm_result=ESC_result
			ypos_track=9
		} else {
			hmm_result=EpiSC_result
			ypos_track=3
		}
		sel=hmm_result[,"Chromosome"]==paste("chr",SE_coordinates[1],sep="") &
		    hmm_result[,"Start"]>=SE_coordinates[2] &
		    hmm_result[,"End"]<=SE_coordinates[3]
	
		meth_SE=hmm_result[sel,"Methylation"]
		state_SE=hmm_result[sel,"State"]
		start_SE=hmm_result[sel,"Start"]
	
		#Methylation track
		scale_track=4
		#baseline of the methylation track
		abline(h=ypos_track,col="darkgrey")
		for (i in 1:length(start_SE)) {
			meth_level=meth_SE[i]
			if (state_SE[i]=="non-methylated") {
				#Inverse methylation level for unmethylated states
				meth_level=-1+meth_level
				meth_color="royalblue2"
			} else {
				meth_color="darkgrey"
			}
			meth_level=ypos_track+scale_track/2*meth_level
			lines(x=c(start_SE[i],start_SE[i]),
			      y=c(ypos_track,meth_level),
			      col=meth_color
			)
		}
	}
	#DM, INT and PU regions
	ypos_track=6
	barheight=0.4
	sel=DM_regions[,1]==paste("chr",SE_coordinates[1],sep="")
	rect(DM_regions[sel,2],ypos_track-barheight/2,DM_regions[sel,3],ypos_track+barheight/2,col="magenta",border="magenta")
	sel=PU_regions[,1]==paste("chr",SE_coordinates[1],sep="")
	rect(PU_regions[sel,2],ypos_track-barheight/2,PU_regions[sel,3],ypos_track+barheight/2,col="green",border="green")
	sel=INT_regions[,1]==paste("chr",SE_coordinates[1],sep="")
	rect(INT_regions[sel,2],ypos_track-barheight/2,INT_regions[sel,3],ypos_track+barheight/2,col="grey",border="grey")

	abline(h=c(1,5,7,11),col="grey")
	axis(side=2,at=c(3,9),labels=c("EpiSC","ESC"),tick=FALSE)
	axis(side=2,
	     at=c(1,3,5,7,9,11),
	     labels=c("-1","0","1","-1","0","1"),
	     las=2,cex.axis=0.8,mgp=c(3, 0.2, 0),tick=FALSE
	)
	dev.off()
}

