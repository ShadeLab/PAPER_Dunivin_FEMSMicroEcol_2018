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
setwd(paste(wd, "/data", sep = ""))

#read in data
names=list.files(pattern="*_matches.txt")
data <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read.delim(X, sep = "\t",header = FALSE))}))

#tidy data so it is readable/ comparable
data.tidy <- data %>%
  separate(V1, into = c("junk", "scaffold"), remove = TRUE, sep = ":") %>%
  select(-junk) %>%
  unique() 

#make a wide dataframe that shows presence of scaffolds in separate genes
data.wide <- dcast(data.tidy, id~scaffold, value.var = "V3")

#replace length with 1 or 0
data.wide[data.wide > 0] <- 1

#remove all columns who only match one contig
number.hits <- colSums(data.wide[2:ncol(data.wide)])
multi.hits <- data.frame(number.hits[number.hits >1])

#remove all scaffolds not present in multiple contigs
data.tidy.matches <- data.tidy[data.tidy$scaffold %in% rownames(multi.hits),]



