#load required packages
library(tidyverse)

#set up environment
setwd("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/diversity_analysis/")
wd <- paste(getwd())

#read in label data
labels <- read_delim("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/phylogenetic_analysis/tree_data/arsM/arsM_0.1_labels_short.txt", delim = ",", col_names = FALSE)
colnames(labels) <- c("Label", "OTU")
labels$OTU <- gsub(" ", "", labels$OTU)

#read in OTU table 
table <- read.delim("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/diversity_analysis/data/0.1_clust/arsM_rformat_dist_0.1.txt", header = TRUE)

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_FULL_map.txt", 
                                    sep=""), sep=" ", header=TRUE))

#read in rplB data
rplB <- read.delim("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/diversity_analysis/output/rplB.summary.scg.txt", header = TRUE, sep = " ")

#add census data for normalization purposes
table.census <- table %>%
  rename(Site = X) %>%
  left_join(rplB, by = "Site") %>%
  select(-Total)

#normalize data
table.normalized <- cbind(table.census$Site, 
                          data.frame(apply(table.census[,2:946], 2, function(x) x/table.census$rplB)))

#rename site column
table.normalized$Site <- table.normalized$`table.census$Site`
table.normalized <- table.normalized[,-1]

#order based on temperature
table.normalized <- table.normalized %>%
  left_join(meta, by = "Site") %>%
  arrange(Site, desc(SoilTemperature_to10cm)) %>%
  select(OTU_001:Site, SoilTemperature_to10cm)


table.normalized$Site <- factor(table.normalized$Site, 
                                levels = table.normalized$Site[order(meta$SoilTemperature_to10cm)])

#transform otu table
table.normalized.t = setNames(data.frame(t(table.normalized[,-c(946,947)])), 
                              table.normalized[,946])
#make OTUs a column
table.normalized.t$OTU <- rownames(table.normalized.t)

#add name to first column
table.normalized.t$OTU <- gsub(" ", "", table.normalized.t$OTU)

#fix OTU names with weird # before them
#table.normalized.t$OTU <- gsub(".OTU_", "OTU_", table.normalized.t$OTU)

#fix label names
#labels$OTU <- gsub("TU_0003", "OTU_0003", labels$OTU)
#labels$OTU <- gsub("TU_0005", "OTU_0005", labels$OTU)
#labels$OTU <- gsub(".OTU_", "OTU_", labels$OTU)

#add leading zero to 4 digits
library(stringr)
rownames(table.normalized.t) <- gsub("OTU_", "", rownames(table.normalized.t))
rownames(table.normalized.t) <- sprintf("%04s", rownames(table.normalized.t))
rownames(table.normalized.t) <- paste("OTU_", rownames(table.normalized.t), sep="")
table.normalized.t$OTU <- rownames(table.normalized.t)
label.abund <- left_join(labels, table.normalized.t, by = "OTU")

#remove OTU information
label.abund <- select(label.abund, -OTU)

#save file
write.csv(label.abund, "/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/phylogenetic_analysis/tree_data/arsM/arsM_abund_label.csv", row.names = FALSE)




