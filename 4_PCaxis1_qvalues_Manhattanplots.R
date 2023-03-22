#############################Plot of PC1 qvalues################################
#Barbara Langille
#March 2023


##Will need bed/fam/bim files 
#I like to use a map file too but you can substitute bim with a few changes to the code
#can get the map/ped files in plink

##this script will calculate qvlaues of PC axis one values - it also picks out outliers
##Manhattan plots


################################################################################




install.packages("devtools")
library("devtools")
install_github("jdstorey/qvalue")


#R packages
library(pcadapt)
library(qvalue)
library(data.table)
library(qqman)
library(ggplot2)


#get into the working directory
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())


#Get the data - add your bed filename here
Datafile <- read.pcadapt(BEDfile,type="bed")

#Rerun with K=1 to get outlier loci associated with axis 1 
pcadapt_1pc <- pcadapt(Datafile, method = "mahalanobis", K=1, min.maf = 0.01)
pvals <- data.frame(pcadapt_1pc$pvalues)

#need to add your map file name here
MAP <- read.delim(MAPfile, stringsAsFactors = F, header = F, sep = "")
PC1_MAP <- data.frame(MAP, pvals)
colnames(PC1_MAP)[1] <- "CHR"
colnames(PC1_MAP)[2] <- "SNP"
colnames(PC1_MAP)[3] <- "X"
colnames(PC1_MAP)[4] <- "POS"
colnames(PC1_MAP)[5] <- "pvalue"

#make the manhattan plot - you will have to change the x and y axis values to suit your data
manhattan(x=PC1_MAP[!is.na(PC1_MAP$CHR),],
          chr="CHR", bp ="POS", p = "pvalue", snp = "SNP",
          genomewideline = F,
          suggestiveline =  F,  #col = c("dodgerblue4", "gray48"),
          xlim=c(16000000,17000000),
          ylim=c(0,60),
          col=c("gray50", "dodgerblue4"),
          highlight=PC1_MAP$SNP[order(PC1_MAP$pvalue, decreasing=F)][1:10]) # change for however many values you want highlighted

#Get outliers based on p-values
alpha <- 0.05
qval <- pcadapt_1pc$pvalues
outliers <- which(qval<alpha)
length(outliers) #number of outliers 


#Get outliers based on q-values
alpha <- 0.05
qval1 <- qvalue(pcadapt_1pc$pvalues)
outliers1 <- which(qval1$qvalues<alpha)
length(outliers1) #number of outliers 


#Create dataframe with qvalues
qvalues <- data.frame(qval1$qvalues) 
colnames(qvalues)[1] <- "qvalue"
#Combine SNP names and qvalue data
results_pc1 <- cbind(PC1_MAP, qvalues) #qvalues are in order of SNPs in map file 

#Percentage of outliers
(nrow(results_pc1[which(results_pc1$qvalue<0.05),])/nrow(results_pc1))*100


#Save data for other analyses
write.table(results_pc1, "PCADAPT_Results_pvalsPC1_K1.txt", 
            col.names = T, row.names = F, quote = F, sep="\t")
write.table(results_pc1[which(results_pc1$qvalue<0.05),], "PCADAPT_outliers_pvalsPC1_K1.txt", 
            col.names = T, row.names = F, quote = F, sep="\t")

#Remove blanks/NAs as won't work in Manhattan plot - you will have to change the x and y axis values
manhattan(x=results_pc1[!is.na(results_pc1$CHR),],
          chr="CHR", bp ="POS", p = "qvalue", snp = "SNP",
          genomewideline = F,
          suggestiveline =  F,  #col = c("dodgerblue4", "gray48"),
          xlim=c(16000000,17000000),
          ylim=c(0,60),col=c("gray50", "dodgerblue4"),
          highlight=results_pc1$SNP[order(results_pc1$qvalue, decreasing=F)][1:10]) # change for however many values you want highlighted


