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

#make organism a character
data$organism <- as.character(data$organism)

#get taxonomy for organisms
data.ncbi <- tax_name(query = data$organism, 
                      get = c("genus", "class", "phylum"), db = "ncbi")

#change query column to "organism"
data.ncbi$organism <- data.ncbi$query

#replace NA in phylum with unknown
data.ncbi$phylum[is.na(data.ncbi$phylum)] = "Unknown"

#call NA class by phyla
data.ncbi$class[is.na(data.ncbi$class)] <- as.character(data.ncbi$phylum[is.na(data.ncbi$class)])

#join taxanomic information with data
data.taxa <- data %>%
  left_join(data.ncbi, by = "organism") %>%
  unique() 


data.phylum <- data.taxa %>%
  group_by(Type, Gene, phylum, Site) %>%
  summarise(Phylum.count = sum(Normalized.Abundance.census))


#prep colors for diversity
color <- c("#FF7F00", "#7570B3", "#CAB2D6", "#FBB4AE", "#F0027F", "#BEBADA", "#E78AC3", "#A6D854", "#B3B3B3", "#386CB0", "#BC80BD", "#FFFFCC", "#BF5B17", "#984EA3", "#CCCCCC", "#FFFF99", "#B15928", "#F781BF", "#FDC086", "#A6CEE3", "#FDB462", "#FED9A6", "#E6AB02", "#E31A1C", "#B2DF8A", "#377EB8", "#FCCDE5", "#80B1D3", "#FFD92F", "#33A02C", "#66C2A5", "#666666", "black", "brown")
(gene.bar.census <- ggplot(data.phylum, aes(x = Site, 
                                            y = Phylum.count, fill = phylum)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = color) +
    theme_classic() +
    facet_wrap(~ Gene, scales = "free_y", ncol = 2))

data.class <- data.taxa %>%
  group_by(Type, Gene, class, Site) %>%
  summarise(Class.count = sum(Normalized.Abundance.census)) %>%
  left_join(data.ncbi, by = "class") %>%
  select(Type:Class.count, phylum) %>%
  unique()

#make class a factor based on phylum
data.class$class <- factor(data.class$class, 
                           data.class$class[order(data.class$phylum)])

#prep colors
n <- 58
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color.class <- print(sample(col_vector, n))


(gene.bar.census.cl <- ggplot(data.class, aes(x = Site, 
                                            y = Class.count, fill = class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color.class) +
    theme_classic() +
    facet_wrap(~ Gene, scales = "free_y", ncol = 2))

(gene.bar.census.cl <- ggplot(data.class, aes(x = Site, 
                                              y = Class.count, color = class)) +
    geom_jitter(stat = "identity", width = 0.05) +
    scale_color_manual(values = color.class) +
    theme_classic() +
    facet_wrap(~ Gene, scales = "free_y", ncol = 2))

data.taxa.classif <- data.taxa %>%
  ungroup() %>%
  group_by(Type, Gene, Site) %>%
  summarise(Abund = sum(Normalized.Abundance.census)) %>%
  left_join(meta, by = "Site")

ggplot(data.taxa.classif, aes(x = Classification, y = Abund)) +
  geom_boxplot() +
  geom_jitter(aes(color = SoilTemperature_to10cm)) +
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                        guide_legend(title="Temperature (°C)")) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2)


ggplot(data.taxa.classif, aes(x = SoilTemperature_to10cm, y=Abund)) +
  geom_smooth(method = "auto") +
  geom_point(aes(shape = Classification)) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2)






