#load required packages
library(tidyverse)

#set up environment
setwd("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/diversity_analysis/")
wd <- paste(getwd())

#read in label data
labels <- read_delim("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/phylogenetic_analysis/aioA/aioA_0.1_labels.txt", delim = ",", col_names = FALSE)
colnames(labels) <- c("Label", "OTU")
labels$OTU <- gsub(" ", "", labels$OTU)

#read in OTU table 
table <- read.delim("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/phylogenetic_analysis/aioA/rformat_dist_0.1.txt", header = TRUE)

#read in census data
census <- read_delim(file = paste(wd, "/data/microbe_census.txt", sep = ""),
                     delim = "\t", col_types = list(col_character(),
                                                    col_number(),
                                                    col_number(), 
                                                    col_number()))
census <- census %>%
  select(Site, GE)

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_JGI_map.txt", 
                                    sep=""), sep=" ", header=TRUE))

#add census data for normalization purposes
table.census <- table %>%
  rename(Site = X) %>%
  left_join(census, by = "Site")

#normalize data
table.normalized <- 100*table.census[2:44]/table.census$GE

#add back site information
table.normalized <- cbind(table.census$Site, table.normalized)

#rename site column
table.normalized$Site <- table.normalized$`table.census$Site`
table.normalized <- table.normalized[,-1]

#order based on temperature
table.normalized <- table.normalized %>%
  left_join(meta, by = "Site") %>%
  arrange(Site, desc(SoilTemperature_to10cm)) %>%
  select(OTU_01:Site, SoilTemperature_to10cm)


table.normalized$Site <- factor(table.normalized$Site, 
                                levels = table.normalized$Site[order(meta$SoilTemperature_to10cm)])

#transform otu table
table.normalized.t = setNames(data.frame(t(table.normalized[,-c(44,45)])), 
                              table.normalized[,44])
#make OTUs a column
table.normalized.t$OTU <- rownames(table.normalized.t)

#add name to first column
table.normalized.t$OTU <- gsub(" ", "", table.normalized.t$OTU)
            
#join abundance with name
label.abund <- left_join(labels, table.normalized.t, by = "OTU")

#remove OTU information
label.abund <- select(label.abund, -OTU)

#save file
write.csv(label.abund, "/Users/dunivint/Documents/ShadeLab/Experiments/Xander/aioA_abund_label.csv", row.names = FALSE)




