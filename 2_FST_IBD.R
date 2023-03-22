######################Fst and Isolation By Distance (IBD)#######################
##Barbara Langille


#will need a vcf file
#This script contains pairwise Fst calculation using stamppFst (StAMPP package)
#This script contains marmap code to get distances between sites
#This script contains Isolation by distance using regression between Fst and distance







#load these packages
library(vcfR)
library(adegenet)
library(StAMPP)
library(dplyr)


#get into the working directory
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())


#start by getting the vcf file in plink v1.9 using the following code
#plink --file <name of ped/map files>  --recode vcf --allow-extra-chr --out <name>

#this is the function - dont need to change anything
Fst_file <- function(VCFfile, FAMfile){
  #read vcf file into into R
  vcf<-read.vcfR(VCFfile, verbose = FALSE)
  #convert vcf into genlight object
  x <- vcfR2genlight(vcf)
  #get the fam file
  FAM <- read.delim(FAMfile, stringsAsFactors = F, header = F, sep = "")
  #add in the individual and population names - needs to be same order as vcf
  x@ind.names <- FAM$V2
  x@pop <- as.factor(FAM$V1)
  #convert the genotype data into allele frequencies
  geno <- stamppConvert(x, type = "genlight")
  #run the pairwise Fst - if your file is too big, run it in arlequin
  FST_table <- stamppFst(geno = geno, nboots = 500, percent = 95, nclusters = 10) 
  #write the values (Fst and p) to tables
  write.csv(FST_table$Fsts, file="FST.csv")
  write.csv(FST_table$Pvalue, file="FST_pvalue.csv")
  
}

#run the function
Fst_file(VCFfile, FAMfile)





################################################################################




##Making a nice Fst figure##




#get into the working directory
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())

#this is the function - dont need to change anything
FST_function <- function(FSTfile){
  #load in data from either arlequin or import the data from above
  #fst_df <- read.delim("EUandNA_Fst.txt", sep = " ", row.names=1) #for arlequin file
  fst_df <- read.csv(FSTfile, header=TRUE, row.names=1) #for data generated above
  ## you will need to change it from a triangle table to long format!!
  ind <- which(lower.tri(fst_df, diag = TRUE), arr.ind = TRUE) ##might need to change between upper.tri and lower.tri depending on what dataset you are using
  nn <- dimnames(fst_df)
  FST <- data.frame(row = nn[[1]][ind[, 1]],
                    col = nn[[2]][ind[, 2]],
                    val = fst_df[ind])
  #write the long format data out
  write.csv(FST, file="FST_longformat.csv")
  #get rid of NAs
  FST[is.na(FST)] <- 0
  #change the significant figs
  FST$val <- formatC(signif(FST$val, digits=2))
  FST$val <- as.numeric(FST$val)
  # Fst italic label for the graph
  fst.label = expression(italic("F")[ST])
  # Extract middle Fst value for gradient argument
  mid = 0.5 # or choose an arbitrary point between 0 and 1
  #make the plot
  ggplot(data = FST, aes(x = row, y = col, fill = val))+
    geom_tile(colour = "black")+
    geom_text(aes(label = val), color="black", size = 3)+
    scale_fill_gradient2(low = "yellow", mid = "pink", high = "red", 
                         midpoint = mid, name = fst.label, 
                         limits = c(0, max(FST$val)), 
                         breaks = c(0.05, 0.25, 0.50, 0.75))+
    scale_x_discrete(expand = c(0,0), position = "bottom")+
    scale_y_discrete(expand = c(0,0), position = "right")+
    theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "right",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 10)
    )
}

#function to run - use the file generated from the step above
FST_function("FST.csv")


################################################################################
################################################################################




# Calculate geographic distances between marine sites using marmap
# Will need this file for IBD below



# Load libraries
library(marmap)
library(Matrix)
library(gdistance)
library(tidyverse)



setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())

#Read in a dataframe with site, latitude, and longitude columns (Lat and Long in header)
#need to make this file 
coords = read.csv("SiteLatLong.csv")

maxlong = max(coords$Long)
minlong = min(coords$Long)
maxlat = max(coords$Lat)
minlat = min(coords$Lat)
# Download bathymetric data
# The higher the resolution the higher accuracy but the slower the run time
# res = 1 is the highest resolution
bathy_data = getNOAA.bathy(lon1 = minlong, lon2 = maxlong, lat1 = minlat, lat2 = maxlat, res = 2, keep = TRUE)
summary(bathy_data)

# Get depth of coordinates
# Must be less than ten metres (<-10) to compute distances (may need to adjust coordinates)
get.depth(bathy_data, x = coords$Long, y = coords$Lat, locator = FALSE)

# Create transition object [long run time]
# Author of marmap recommendations:
# Use a minimum depth of -10 to avoid path crossing land masses
# Use a maximum depth of -200 to limit paths to continental shelf
trans1 = trans.mat(bathy_data, min.depth = -10, max.depth = NULL)

# Compute least-cost paths [long run time]
lc_paths = lc.dist(trans1, subset(coords, select = c("Long","Lat")), res = c("path"))

# Compute least-cost distances (km) matrix
lc_dist = lc.dist(trans1, subset(coords, select = c("Long","Lat")), res = "dist")

# Convert to matrix, rename columns and rows, and export as csv file
lc_mat = as.matrix(lc_dist)
colnames(lc_mat) = as.vector(coords$Pop)
rownames(lc_mat) = as.vector(coords$Pop)
lc_mat
write.csv(lc_mat, file="marmap_distances_km.csv")

# Plot paths on a map
# Visually check that no path overlaps land
map <- plot.bathy(bathy_data, image= TRUE, land = TRUE, n = 0,
                  bpal = list(c(0, max(bathy_data), "grey"),
                              c(min(bathy_data), 0, "royalblue")))
map + lapply(lc_paths, lines, col = "orange", lwd = 2, lty = 1)




################################################################################
################################################################################




####Isolation by distance (IBD)####




#Load these packages
library(geosphere)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(Matrix)
library(viridis)

#get into the working directory
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())

#Need to start by calculating the geographic distance between populations
#DO the second function below if you are using marmap distances (from above)
##import csv with lat/lon - make sure same order as pops in fst file
IBD_crowflies <- function(METADATAfile, FSTfile_long, FAMfile){
  geodata <- read.csv(METADATAfile, header=TRUE, row.names = 1)
  geodata1 <- as.matrix(geodata)
  #calculate the distance using an as the crow flies measure
  #you can get a more accurate distance using marmap but this will do for non-cyclical areas
  dist <- distm(geodata1, fun = distHaversine)
  geodata_names <- read.csv(METADATAfile, header=TRUE)
  rownames(dist) <- c(geodata_names$Pop)
  colnames(dist) <- c(geodata_names$Pop)
  ##transform dist file into long format
  ind <- which(lower.tri(dist, diag = TRUE), arr.ind = TRUE)
  nn <- dimnames(dist)
  new_dist <- data.frame(row = nn[[1]][ind[, 1]],
                         col = nn[[2]][ind[, 2]],
                         val = dist[ind])
  #write the long format for the dist file 
  write.csv(new_dist,"dist_longformat.csv") 
  colnames(new_dist)[3] <- "distance"
  #get distance in kms
  new_dist1 <- data.frame((new_dist$distance)/1000)
  ##Bring the FST data back in
  FSTfile <- read.csv(FSTfile_long, header=TRUE, row.names = 1)
  FSTfile[is.na(FSTfile)] <- 0
  #combine FST and Geo file
  colnames(FSTfile)[3] <- "FST" 
  FST_GEO <- data.frame(c(FSTfile, new_dist1))
  colnames(FST_GEO)[4] <- "distance"
  #add in a proper regression line
  fit1 <- lm(FST_GEO$FST ~ FST_GEO$distance)
  print(summary(fit1)) #view the output and get R2 and adjusted R2
  #get the fam file
  FAM <- read.delim(FAMfile, stringsAsFactors = F, header = F, sep = "")
  #add in the population information
  fit1$Pop <- as.factor(FAM$V1)
  ##plot regressing FST and geographic distance
  ggplotRegression <- function (fit1) {
    require(ggplot2)
    ggplot(fit1$model, aes(x = fit1$model$`FST_GEO$dist`, 
                           y = fit1$model$`FST_GEO$FST`)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "red", se=FALSE) +
      #scale_x_continuous(limits=c(0,10400), expand = c(0, 0)) + ##need to change distance based on what you have
      #scale_y_continuous(limits=c(0,0.18), expand = c(0, 0)) + ##need to change fst based on what you have
      theme_classic() +
      labs(title = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5)), 
           x ="Geographical Distance (km)", y ="FST")
  }
  #get the plot
  ggplotRegression(lm(FST_GEO$FST ~ FST_GEO$dist))
  
}

#run the function for an as the crow flies measure - Fst file was generated above
IBD_crowflies(METADATAfile, "FST_longformat.csv", FAMfile)

##dont forget to save your plot
##summary stats should be displayed in console


################################################################################


#use this function if you are using the marmap distances
IBD_marmap <- function(DISTfile, FSTfile_long, FAMfile){
  dist <- read.csv(DISTfile, header=TRUE, row.names=1) #for data generated above
  ## you will need to change it from a triangle table to long format!!
  ind <- which(lower.tri(dist, diag = TRUE), arr.ind = TRUE) ##might need to change between upper.tri and lower.tri depending on what dataset you are using
  nn <- dimnames(dist)
  dist_longformat <- data.frame(row = nn[[1]][ind[, 1]],
                                col = nn[[2]][ind[, 2]],
                                val = dist[ind])
  #write the long format for the dist file 
  write.csv(dist_longformat,"dist_marmap_longformat.csv") 
  colnames(dist_longformat)[3] <- "distance"
  new_dist1 <- data.frame(dist_longformat$distance)
  ##Bring the FST data back in
  FSTfile <- read.csv(FSTfile_long, header=TRUE, row.names = 1)
  FSTfile[is.na(FSTfile)] <- 0
  #combine FST and Geo file
  colnames(FSTfile)[3] <- "FST" 
  FST_GEO <- data.frame(c(FSTfile, new_dist1))
  colnames(FST_GEO)[4] <- "distance"
  #add in a proper regression line
  fit1 <- lm(FST_GEO$FST ~ FST_GEO$distance)
  print(summary(fit1)) #view the output and get R2 and adjusted R2
  #get the fam file
  FAM <- read.delim(FAMfile, stringsAsFactors = F, header = F, sep = "")
  #add in the population information
  fit1$Pop <- as.factor(FAM$V1)
  ##plot regressing FST and geographic distance
  ggplotRegression <- function (fit1) {
    require(ggplot2)
    ggplot(fit1$model, aes(x = fit1$model$`FST_GEO$dist`, 
                           y = fit1$model$`FST_GEO$FST`)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "red", se=FALSE) +
      #scale_x_continuous(limits=c(0,10400), expand = c(0, 0)) + ##need to change distance based on what you have
      #scale_y_continuous(limits=c(0,0.18), expand = c(0, 0)) + ##need to change fst based on what you have
      theme_classic() +
      labs(title = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5)), 
           x ="Geographical Distance (km)", y ="FST")
  }
  #get the plot
  ggplotRegression(lm(FST_GEO$FST ~ FST_GEO$dist))
  
}

#run the function - use the fst file generated above
IBD_marmap(DISTfile, "FST_longformat.csv", FAMfile)


##dont forget to save your plot
##summary stats should be displayed in console

################################################################################
