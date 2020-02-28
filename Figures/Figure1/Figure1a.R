#setwd("C:/Users/LucJouneau/Documents/INRA/Alice/Article Azuara/article/Figure1")

#-----------------
# Methylation call
#-----------------

SE_coordinates=c(4,55469259,55491081)
chr=SE_coordinates[1]
maxGap=999

tab=read.table(file=paste("data/cnv384_CpG_methcounts_chr",chr,".bed",sep=""),sep="\t",header=F)

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

#-----------------------
# Transcription factors
#-----------------------
tfs=list()
for (tf_name in c("OCT4","SOX2","NANOG","ESRRB","NCOA3","KLF4","MED1","P300","SMAD3","STAT3","H3K27ac","H3K4me1")) {
	tf_file=list.files(path="bed_cistrome_mm9",pattern=paste(tf_name,".*",sep=""),full.names=TRUE)
	tfs[[tf_name]]=read.table(file=tf_file,header=FALSE,sep="\t")
}

#----------------------------------
# Restrict data to KLF4-SE region
#----------------------------------
sel_KLF4=hmm_result[,"Chromosome"]==paste("chr",SE_coordinates[1],sep="") &
	 hmm_result[,"Start"]>=SE_coordinates[2] &
	 hmm_result[,"End"]<=SE_coordinates[3]

meth_KLF4=hmm_result[sel_KLF4,"Methylation"]
state_KLF4=hmm_result[sel_KLF4,"State"]
start_KLF4=hmm_result[sel_KLF4,"Start"]

#-----------------
# Figure 1
#-----------------
jpeg("Figure1a.jpg")
mar=par("mar")
mar_svg=mar
mar[1]=0.5
mar[3]=6
par("mar"=mar)
plot(x=c(),y=c(),
     xlim=c(SE_coordinates[2],SE_coordinates[3]),ylim=c(0,15),
     main="Klf4-associated super-enhancer",
     xlab="",ylab="",xaxt="n",yaxt="n"
)

#Methylation track
scale_track=2; ypos_track=14
#baseline of the methylation track
abline(h=ypos_track,col="darkgrey")
for (i in 1:length(start_KLF4)) {
	meth_level=meth_KLF4[i]
	if (state_KLF4[i]=="non-methylated") {
		#Inverse methylation level for unmethylated states
		meth_level=-1+meth_level
		meth_color="royalblue2"
	} else {
		meth_color="darkgrey"
	}
	meth_level=ypos_track+scale_track/2*meth_level
	lines(x=c(start_KLF4[i],start_KLF4[i]),
	      y=c(ypos_track,meth_level),
	      col=meth_color
	)
}

abline(h=12.5,col="grey")
scale_track=1;ypos_track=12
for (tf_name in names(tfs)) {
	tf=tfs[[tf_name]]
	sel_tf= tf[,1]==paste("chr",SE_coordinates[1],sep="") &
		tf[,2]>=SE_coordinates[2] &
		tf[,3]<=SE_coordinates[3]
	if (tf_name %in% c("P300","H3K27ac","H3K4me1")) {
		tf_color="green"
		bar_height=0.3
		line_width=1
	} else {
		tf_color="red"
		bar_height=0.2
		line_width=2
	}
	rect(tf[sel_tf,2],ypos_track-bar_height/2,
	     tf[sel_tf,3],ypos_track+bar_height,
	     border=tf_color,col=tf_color,lwd=2
	)
	ypos_track=ypos_track-1
}
axis(side=3,cex.axis=0.8,tcl=-0.25)
axis(side=2,
     at=c(1:length(tfs),13:15),
     labels=c(rev(names(tfs)),"-1","0","1"),
     font=rep(0,length(tfs),0,2,0),
     las=1,cex.axis=0.8,tick=FALSE,
     mgp=c(3, 0.2, 0)
)
axis(side=2,at=14,font=2,labels="mCpG   ",las=1,cex.axis=0.9,tick=FALSE)

par("mar"=mar_svg)
dev.off()
