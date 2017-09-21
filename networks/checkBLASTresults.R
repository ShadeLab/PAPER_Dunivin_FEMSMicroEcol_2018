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
setwd(paste(wd, "/../networks/blast_data", sep = ""))

#read in data
names=list.files(pattern="*_blast.txt")
data <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read.delim(X, sep = "\t",header = FALSE))}))

#tidy data
data$id <- gsub("arsC_", "arsC", data$id)
data.tidy <- data %>%
  separate(id, into = c("Gene", "Clust", "Site"), sep = "_", remove = FALSE)

#make a wide dataframe that shows presence of scaffolds in separate genes
data.wide <- dcast(data.tidy, id~V1+Site, value.var = "V3")

#replace length with 1 or 0
data.wide[data.wide > 0,c(2:ncol(data.wide))] <- 1

#remove all columns who only match one contig
number.hits <- colSums(data.wide[2:ncol(data.wide)])
multi.hits <- data.frame(number.hits[number.hits >1])
multi.hits$scaffold <- rownames(multi.hits)
multi.hits <- multi.hits %>%
  separate(scaffold, into = c("scaffold", "Site"), sep = "_")
  

#remove all scaffolds not present in multiple contigs
data.tidy.matches <- data.tidy[data.tidy$V1 %in% multi.hits$scaffold,]



