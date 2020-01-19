######################################################################
# Wout Megchelenbrink
# Spearman rank correlation between all capture Hi-C experiments
#
#  Jan 17, 2020
#
######################################################################
rm(list=ls())
source("include/style.R")
source("include/include.R")

# Read the transcript-SElet C-HiC data for the Dpn2 and Hind3 based experiments
dpn2  <- fread("DATA/CHIC_tx_SElet_Dpn2.tsv")
hind3 <- fread("DATA/CHIC_tx_SElet_Hind3.tsv")

# Aggregate on the superenhancer level
dpn2 <- dpn2[, lapply(.SD, sum), by=.(gene_name, transcript_id, se_id, nb_pmt_joshi, nb_pmt_sahlen), .SDcols = c("reads_2i_1", "reads_2i_2", "reads_ser_1", "reads_ser_2", "reads_sahlen")]
hind3 <- hind3[, lapply(.SD, sum), by=.(gene_name, transcript_id, se_id, nb_pmt_novo), .SDcols = c("reads_novo_ser")]

# Merge Dpn2 and Hind3 data sets
reads <- merge(dpn2, hind3, by=c("gene_name", "transcript_id", "se_id"))

# This seems like quite a tedious way to get the correlation matrix, but since it ain't broken, let's not fix it ;-)
cor.se <- matrix(NA, nrow=6, ncol=6)
cor.se[1,1] <- reads[nb_pmt_joshi > 0, cor(reads_2i_1, reads_2i_1, method = "spearman")]          
cor.se[2,1] <- reads[nb_pmt_joshi > 0, cor(reads_2i_1, reads_2i_2, method = "spearman")]          
cor.se[3,1] <- reads[nb_pmt_joshi > 0, cor(reads_2i_1, reads_ser_1, method = "spearman")]          
cor.se[4,1] <- reads[nb_pmt_joshi > 0, cor(reads_2i_1, reads_ser_2, method = "spearman")]          
cor.se[5,1] <- reads[nb_pmt_joshi > 0 & nb_pmt_sahlen > 0, cor(reads_2i_1, reads_sahlen, method = "spearman")] 
cor.se[6,1] <- reads[nb_pmt_joshi > 0 & nb_pmt_novo > 0, cor(reads_2i_1, reads_novo_ser, method = "spearman")] 

cor.se[1,2] <- reads[nb_pmt_joshi > 0, cor(reads_2i_2, reads_2i_1, method = "spearman")]          
cor.se[2,2] <- reads[nb_pmt_joshi > 0, cor(reads_2i_2, reads_2i_2, method = "spearman")]          
cor.se[3,2] <- reads[nb_pmt_joshi > 0, cor(reads_2i_2, reads_ser_1, method = "spearman")]          
cor.se[4,2] <- reads[nb_pmt_joshi > 0, cor(reads_2i_2, reads_ser_2, method = "spearman")]          
cor.se[5,2] <- reads[nb_pmt_joshi > 0 & nb_pmt_sahlen > 0, cor(reads_2i_2, reads_sahlen, method = "spearman")] 
cor.se[6,2] <- reads[nb_pmt_joshi > 0 & nb_pmt_novo > 0, cor(reads_2i_2, reads_novo_ser, method = "spearman")] 

cor.se[1,3] <- reads[nb_pmt_joshi > 0, cor(reads_ser_1, reads_2i_1, method = "spearman")]          
cor.se[2,3] <- reads[nb_pmt_joshi > 0, cor(reads_ser_1, reads_2i_2, method = "spearman")]          
cor.se[3,3] <- reads[nb_pmt_joshi > 0, cor(reads_ser_1, reads_ser_1, method = "spearman")]          
cor.se[4,3] <- reads[nb_pmt_joshi > 0, cor(reads_ser_1, reads_ser_2, method = "spearman")]          
cor.se[5,3] <- reads[nb_pmt_joshi > 0 & nb_pmt_sahlen > 0, cor(reads_ser_1, reads_sahlen, method = "spearman")] 
cor.se[6,3] <- reads[nb_pmt_joshi > 0 & nb_pmt_novo > 0, cor(reads_ser_1, reads_novo_ser, method = "spearman")] 

cor.se[1,4] <- reads[nb_pmt_joshi > 0, cor(reads_ser_2, reads_2i_1, method = "spearman")]          
cor.se[2,4] <- reads[nb_pmt_joshi > 0, cor(reads_ser_2, reads_2i_2, method = "spearman")]          
cor.se[3,4] <- reads[nb_pmt_joshi > 0, cor(reads_ser_2, reads_ser_1, method = "spearman")]          
cor.se[4,4] <- reads[nb_pmt_joshi > 0, cor(reads_ser_2, reads_ser_2, method = "spearman")]          
cor.se[5,4] <- reads[nb_pmt_joshi > 0 & nb_pmt_sahlen > 0, cor(reads_ser_2, reads_sahlen, method = "spearman")] 
cor.se[6,4] <- reads[nb_pmt_joshi > 0 & nb_pmt_novo > 0, cor(reads_ser_2, reads_novo_ser, method = "spearman")] 

cor.se[1,5] <- reads[nb_pmt_sahlen > 0 & nb_pmt_joshi > 0, cor(reads_sahlen, reads_2i_1, method = "spearman")] 
cor.se[2,5] <- reads[nb_pmt_sahlen > 0 & nb_pmt_joshi > 0, cor(reads_sahlen, reads_2i_2, method = "spearman")] 
cor.se[3,5] <- reads[nb_pmt_sahlen > 0 & nb_pmt_joshi > 0, cor(reads_sahlen, reads_ser_1, method = "spearman")] 
cor.se[4,5] <- reads[nb_pmt_sahlen > 0 & nb_pmt_joshi > 0, cor(reads_sahlen, reads_ser_2, method = "spearman")] 
cor.se[5,5] <- reads[nb_pmt_sahlen > 0 , cor(reads_sahlen, reads_sahlen, method = "spearman")] 
cor.se[6,5] <- reads[nb_pmt_sahlen > 0 & nb_pmt_novo > 0, cor(reads_sahlen, reads_novo_ser, method = "spearman")] 

cor.se[1,6] <- reads[nb_pmt_novo > 0 & nb_pmt_joshi > 0, cor(reads_novo_ser, reads_2i_1, method = "spearman")] 
cor.se[2,6] <- reads[nb_pmt_novo > 0 & nb_pmt_joshi > 0, cor(reads_novo_ser, reads_2i_2, method = "spearman")] 
cor.se[3,6] <- reads[nb_pmt_novo > 0 & nb_pmt_joshi > 0, cor(reads_novo_ser, reads_ser_1, method = "spearman")] 
cor.se[4,6] <- reads[nb_pmt_novo > 0 & nb_pmt_joshi > 0, cor(reads_novo_ser, reads_ser_2, method = "spearman")] 
cor.se[5,6] <- reads[nb_pmt_novo > 0 & nb_pmt_sahlen > 0, cor(reads_novo_ser, reads_sahlen, method = "spearman")] 
cor.se[6,6] <- reads[nb_pmt_novo > 0, cor(reads_novo_ser, reads_novo_ser, method = "spearman")] 

# Round and annotate rows/cols
cor.se <- data.matrix(round(cor.se,2))
rownames(cor.se) <- c("2i 1", "2i 2", "ser 1", "ser 2", "sahlen", "novo")
colnames(cor.se) <- c("2i 1", "2i 2", "ser 1", "ser 2", "sahlen", "novo")

# Make the correlation plot (correlation.plot is defined in "include/include.R")
pdf("IMG/FigS3A.pdf", width = 5.5, height = 5)
correlation.plot(cor.se, title="Rs")
dev.off()
