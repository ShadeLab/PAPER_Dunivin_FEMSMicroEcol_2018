##########################
#READ IN DATA, SET UP ENV#
##########################

#read dependencies
library(phyloseq)
library(vegan)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(taxize)


#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis
wd <- print(getwd())

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_JGI_map.txt", 
                                    sep=""), sep=" ", header=TRUE))

####################################
#EXAMINE AsRG across chronosequence#
####################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/output", sep = ""))

#read in abundance data
names <- list.files(pattern="*summary.txt")
data <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read.table(X, sep = " " , header = TRUE))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#summarise data to get number of genes per gene per site
summarised <- data %>%
  group_by(Type, Gene, Site) %>%
  summarise(Count = sum(Abundance), 
            Count.rplB = sum(Normalized.Abundance.rplB),
            Count.census = sum(Normalized.Abundance.census))

#order based on temperature
summarised$Site <- factor(summarised$Site, 
                          summarised$Site[order(meta$SoilTemperature_to10cm)])

#plot gene data
(gene.plot.census <- ggplot(summarised, aes(x = Site, y = Count.census)) +
    geom_bar(stat = "identity", aes(fill = Gene)) +
    ylab("Gene count (normalized to genome equivalents)") +
    scale_fill_brewer(type = "qualitative", palette = "Paired") +
    theme_classic())

(gene.plot.rplB <- ggplot(summarised, aes(x = Site, y = Count.rplB)) +
    geom_bar(stat = "identity", aes(fill = Gene)) +
    ylab("Gene count (normalized to genome equivalents)") +
    scale_fill_brewer(type = "qualitative", palette = "Paired") +
    theme_classic())

#group data
data <-  data %>%
  group_by(Type, Gene, Site)

#order based on temperature
data$Site <- factor(data$Site, 
                    data$Site[order(data$Temp)])

(gene.bar.censsu <- ggplot(data, aes(x = Site, y = Normalized.Abundance.census)) +
  geom_jitter(aes(color = organism))) +
  facet_wrap(~ Gene, scales = "free_y")












