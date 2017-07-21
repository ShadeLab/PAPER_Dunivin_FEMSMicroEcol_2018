#load dependencies 
library(phyloseq)
library(vegan)
library(tidyverse)
library(reshape2)
library(RColorBrewer)

#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis
wd <- print(getwd())

#setwd to diversity analysis
setwd(paste(wd, "/networks", sep = ""))

#now change wd obj
wd <- print(getwd())

#read in data
names=list.files(pattern="*_blast.txt")
data <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read.delim(X, sep = "\t",header = FALSE))}))

#check how many are unique
length(unique(data$V2))
#==315/317; means 1 duplicate sequence

#Examine which scafold is not unique
duplicated(data$V2)
#==scaffold00202; both arsM; both low %id hits