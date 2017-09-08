#load required packages
library(tidyverse)

#set up environment
wd <- paste(getwd())

#read in OTU table 
table <- read.delim(paste(wd, "/diversity_analysis/data/tetW_rformat_dist_0.1.txt", sep = ""), header = TRUE)
table$X <- gsub("cen", "Cen", table$X)

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/diversity_analysis/data/Centralia_FULL_map.txt", sep=""), sep=" ", header=TRUE))

#read in rplB data
rplB <- read.delim(paste(wd, "/diversity_analysis/output/rplB.summary.scg.txt", sep = ""), header = TRUE, sep = " ")

#add census data for normalization purposes
table.census <- table %>%
  rename(Site = X) %>%
  left_join(rplB, by = "Site") %>%
  rename(rplB = rowSums.rplB.)

#normalize data
table.normalized <- cbind(Site = table.census$Site, 
                          data.frame(apply(table.census[,2:ncol(table.census)], 2, function(x) x/table.census$rplB)))

#transform otu table
table.normalized.t = setNames(data.frame(t(table.normalized[,-c(1,ncol(table.normalized))])), 
                              table.normalized[,1])

#add leading zero to 4 digits
library(stringr)
rownames(table.normalized.t) <- gsub("OTU_", "", rownames(table.normalized.t))
rownames(table.normalized.t) <- sprintf("%04s", rownames(table.normalized.t))
rownames(table.normalized.t) <- paste("OTU_", rownames(table.normalized.t), sep="")

#save file
write.csv(table.normalized.t, paste(wd, "/phylogenetic_analysis/tree_data/tetW_abund_label.csv", sep = ""),row.names = TRUE, quote = FALSE)

print(max(rowSums(table.normalized.t)))
print(min(rowSums(table.normalized.t)))
print(mean(rowSums(table.normalized.t)))



