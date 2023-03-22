################################################################################
################################################################################


####admixture and structure####
##### hybriddetective #####
####     pie charts    ####
#######Barbara Langille########


#below code is pulled from Bradbury lab members and myself


#This script contains admixture plots using qmatrix data (ran admixture in ComputeCanada first)
#####starts at line: 32
#This script contains structure plots using snmf
#####starts at line: 73
#this script contains hybriddetective which was used to simulate aquaculture samples for admixture
#####starts at line: 224
#this script contains code for making pie charts from admixture proportions for populations
#####starts at line: 311



################################################################################
################################################################################





####admixture####




#get into the working directory
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())

#run admixture (I use ComputeCanada) and get a .Q file for the correct clusters

#load in the Q-matrix of the data and write it to a csv if not in the order that 
#you want to use it as - if its correct then forget this step
q_matrix <- read.table("FILEname.Q")
write.csv(q_matrix, "q_matrix.csv")

#put in order that you want and re-load it in - can only have the number of 
#columns associated with the clusters - no individual information, so be aware 
#of the order you are choosing
q_matrix <- read.csv("q_matrix.csv")

#if you want it in order based on admixture proportions use this code
#add in more columns depending on how many clusters you have
ord <- q_matrix[order(q_matrix$V1, q_matrix$V2),] 

#order by populations that I chose (vaguely by latitude) - change col for how 
#many clusters there are
admixture_plot <- barplot(t(as.matrix(q_matrix)), space=c(0.2), col=rainbow(4), xlab= "Individual", ylab="Ancestry", border=NA, las=2)




################################################################################
################################################################################





####structure####




#Load in the packages
library(LEA)
library(tidyverse)
library(reshape2)
library(viridis)

##Load in these for later
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

#get into the working directory
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())

#load in a structure file - made in PGDSpider 
#format = 2 means there is 2 rows of data per each individual
#if you already have the .geno file then go to the next stage
struct2geno(file = "FILEname.str", TESS = TRUE, diploid = TRUE, FORMAT = 2,
            extra.row = 0, extra.col = 2, output = "FILEname.geno")

#figure out the correct number of clusters - may need to change K value
obj.snmf = snmf("FILEname.geno", K = 10, ploidy = 2, entropy = T,
                alpha = 100, project = "new")

#look at all the values of K to pick the best one
plot(obj.snmf, col = "blue", cex = 1.4, pch = 19)

##get qmatrix of the best K value
qmatrix = Q(obj.snmf, K = 2)

#need to format my qvalues
apply(qmatrix, 1, which.max) %>%
  as_tibble() %>%
  dplyr::rename(Genetic_group = value) %>%
  write_tsv('FILEname_group.txt')

#this is the long version - there is a smaller one below. only run one
#change it for how many groups you have
as_tibble(qmatrix) %>%
  dplyr::rename(Q1 = V1,
                Q2 = V2,
                Q3 = V3,
                Q4 = V4,
                Q5 = V5,
                Q6 = V6,
                Q7 = V7,
                Q8 = V8,
                Q9 = V9,
                Q10 = V10,
                Q11 = V11,
                Q12 = V12,
                Q13 = V13,
                Q14 = V14,
                Q15 = V15,
                Q16 = V16) %>%
  write_csv('FILEname_qvalues_k16.csv')

#this is the same as above but its easier to have one thats smaller instead of 
#altering the upper one
as_tibble(qmatrix) %>%
  dplyr::rename(Q1 = V1,
                Q2 = V2) %>%
  write_csv('FILEname_qvalues_k2.csv')

##read in the group data generated above
group = read_tsv('FILEname_group.txt')

#load in the extra variables
environment_data = read_csv('METADATAfile.csv') %>% 
  dplyr::rename(FID = FamID, Lat = Latitude, Long = Longitude, ID = SampleID) %>% 
  dplyr:: select(FID,
                 ID,
                 Lat, 
                 Long)

#read in the qvalues
qvalues = read_csv('FILEname_qvalues_k2.csv')

#put it all together
snmf_data = bind_cols(environment_data, group, 
                      qvalues) %>% 
  arrange(Lat) %>% 
  na.omit()

#write out the structure file
write.csv(snmf_data, "FILEname_STRUCTURE.csv")

snmf_melted = melt(snmf_data, 
                   id.vars = c('FID', 
                               'ID', 
                               'Genetic_group', 
                               'Lat',
                               'Long')) %>% 
  as_tibble()

## ggplot of ancestry proportion graph from snmf!!!!
theme_set(theme_bw())


## snmf latitude plot
ggplot(data = snmf_melted, 
                       aes(x = reorder(ID, Lat),
                           y = value, 
                           fill = variable, 
                           group = Lat))+
  geom_bar(stat = "identity", 
           width = 1)+
  scale_fill_manual(values = magma(n = 2))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 6,
        #                            color = 'black'),
        legend.position = 'none')+
  #scale_x_discrete(breaks = breaks, labels = labels)+ #I cant get this bit with labels and breaks to work
  scale_y_continuous(expand = c(0,0))



#this will show you the order they are in
unique(snmf_melted$FID)


#look at this value in a barplot
barplot(t(snmf_data$...5), col = c("blue","brown","lightgreen","red", "purple", 
                                   "yellow", "green", "orange", "pink", "white", 
                                   "lightblue", "black"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")





################################################################################





####hybriddetective
### how to center data on an ancestor (in this case aquaculture fish)





library(devtools)
devtools::install_github("rystanley/genepopedit") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load(parallel, plyr, stringr, tidyr)
devtools::install_github("bwringe/parallelnewhybrid")

devtools::install_github("bwringe/hybriddetective")

library(genepopedit)
library(parallelnewhybrid)
library(hybriddetective)
library(adegenet)

#get into the working directory
setwd("~/DIRECTORY")

#run hybriddetective - need an older version of R for this to work (4.0.2 works)
freqbasedsim_AlleleSample(GPD = "GENEPOPfile.txt", pop.groups = c("AQ", "NL"))



################################################################################




#get qvalues in snmf



#Load in the packages
library(LEA)
library(tidyverse)
library(reshape2)
library(viridis)

##Load in these for later
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

#get into the working directory
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())

#load in a structure file - made in PGDSpider 
#format = 2 means there is 2 rows of data per each individual
#if you already got the .geno file then go to the next stage
struct2geno(file = "FILEname.str", TESS = TRUE, diploid = TRUE, FORMAT = 2,
            extra.row = 0, extra.col = 2, output = "FILEname.geno")

#figure out the correct number of clusters from 1 to 20
obj.snmf = snmf("FILEname.geno", K = 1:7, ploidy = 2, entropy = T,
                alpha = 100, project = "new")

#look at all the values of K to pick the best one
plot(obj.snmf, col = "blue", cex = 1.4, pch = 19)

##get qmatrix of the best K value
qmatrix = Q(obj.snmf, K = 4)

#many clusters there are
admixture_plot <- barplot(t(as.matrix(qmatrix)), space=c(0.2), col=rainbow(4), 
                          xlab= "Individual", ylab="Ancestry", border=NA, las=2)

#get file of admixture proportions - need this for pie charts below - change Q if you need
as_tibble(qmatrix) %>%
  dplyr::rename(Q1 = V1,
                Q2 = V2) %>%
  write_csv('FILEname_qvalues_k2.csv')



################################################################################




#easy pie charts



##Make a pie chart with admixture proportions

library(ggplot2)

#read in the data - need to sum up the amount of each admixture proportion for each site
Pie <- read.csv("FILEname_piecharts.csv", header = TRUE)

ggplot(Pie, aes(x = "", y = CNR3, fill = Qvalue)) +
  geom_col() +
  coord_polar(theta = "y")



