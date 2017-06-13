#load packages
library(dplyr)
library(reshape2)
library(ggplot2)

#setwd
wd <- paste(getwd())

##################
#set up COGs data#
##################

#read in data
cogs <- read.delim(file = paste(wd, "/COGs/data/COG_data.txt", sep=""))

#call rplB 
rplB <- "COG0090"

#select rplB
rplB <- cogs[which(cogs$COGID %in% rplB),]

#tidy rplB data
rplB <- rplB %>%
  melt(id.vars = c("COGID", "Func_name"), variable.name = "Site", value.name = "Count") %>%
  rename(rplB = Count)

#tidy COG data
cogs <- melt(cogs, id.vars = c("COGID", "Func_name"), variable.name = "Site", value.name = "Count")

#join with rplB information
cogs.rplB <- cogs %>%
  left_join(rplB, by = "Site") %>%
  mutate(Normalized.Abundance.rplB = Count/rplB, 
         Method = "COGs") %>%
  rename(Gene = Func_name.x, Abundance = Count) %>%
  select(Gene, Site, Abundance, Normalized.Abundance.rplB, Method)

#change site from "C" to "Cen" so it matches metadata
cogs.rplB$Site <- gsub("C", "Cen", cogs.rplB$Site)

#list AsRG of interest
asrg <- c("Arsenite efflux pump ArsB, ACR3 family", "Arsenate reductase and related proteins, glutaredoxin family")

#Select only AsRG of interest
cogs.asrg <- cogs.rplB[which(cogs.rplB$Gene %in% asrg),]

#change gene names to match xander
cogs.asrg$Gene <- gsub("Arsenite efflux pump ArsB, ACR3 family", "acr3",
                       cogs.asrg$Gene)

cogs.asrg$Gene <- gsub("Arsenate reductase and related proteins, glutaredoxin family", 
                       "arsCglut",
                       cogs.asrg$Gene)

#####################
#read in xander data#
#####################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/diversity_analysis/data", sep = ""))

#read in abundance data
names <- list.files(pattern="*_45_taxonabund.txt")
xander <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read_delim(X, delim = "\t"))}))

#move back up a directory to proceed with analysis
setwd("../../")
wd <- print(getwd())

#split columns and tidy dataset
xander <- xander %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = junk, into = c("Gene", "junk"), sep = "_45_", remove = TRUE)

#remove awkward _ in Gene column
xander$Gene <- gsub("_", "", xander$Gene)

#change site from "cen" to "Cen" so it matches metadata
xander$Site <- gsub("cen", "Cen", xander$Site)

#separage out rplB data (not needed for gene-centric analysis)
rplB.xander <- xander[which(xander$Gene == "rplB"),]
xander <- xander[-which(xander$Gene == "rplB"),]

#split columns 
rplB.xander <- rplB.xander %>%
  select(Site, Taxon:Fraction.Abundance) %>%
  group_by(Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
rplB.xander$Fraction.Abundance <- as.numeric(rplB.xander$Fraction.Abundance)
rplB.xander$Abundance <- as.numeric(rplB.xander$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.rplB <- rplB.xander %>%
  summarise(Total = sum(Fraction.Abundance), rplB = sum(Abundance))

#Tidy gene data
xander.tidy <- xander %>%
  separate(col = Taxon, into = c("Code", "Organism"), sep = "organism=") %>%
  separate(col = Organism, into = c("Organism", "Definition"), sep = ",definition=") %>%
  select(Site, Gene, Organism:Fraction.Abundance) %>%
  group_by(Gene, Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
xander.tidy$Fraction.Abundance <- as.numeric(xander.tidy$Fraction.Abundance)
xander.tidy$Abundance <- as.numeric(xander.tidy$Abundance)


#make column for organism name and join with microbe census data and normalize to it
xander.abundance <- xander.tidy %>%
  left_join(summarised.rplB, by = "Site") %>%
  mutate( Normalized.Abundance.rplB = Abundance / rplB, 
          Method = "Xander") %>%
  select(Gene, Site, Abundance, Normalized.Abundance.rplB, Method) %>%
  group_by(Method, Gene, Site) %>%
  summarise(Abundance = sum(Abundance))
  
xander.N.abundance <- xander.tidy %>%
  left_join(summarised.rplB, by = "Site") %>%
  mutate( Normalized.Abundance.rplB = Abundance / rplB, 
          Method = "Xander") %>%
  select(Gene, Site, Abundance, Normalized.Abundance.rplB, Method) %>%
  group_by(Gene, Site) %>%
  summarise(Normalized.Abundance.rplB = sum(Normalized.Abundance.rplB)) %>%
  select(Gene, Site, Normalized.Abundance.rplB)

#join together summarised data
xander.annotated <-  xander.abundance %>%
  left_join(xander.N.abundance, by = c("Gene", "Site")) 

#list asrg of interest
asrg.xander <- c("acr3", "arsCglut")

#select xander asrg of interest
xander.asrg <- xander.annotated[which(xander.annotated$Gene %in% asrg.xander),]

#join together datasets
xander.asrg <- data.frame(xander.asrg)
cogs.asrg <- data.frame(cogs.asrg)
mixed <- rbind(xander.asrg, cogs.asrg)

#plot data
ggplot(mixed, aes(x = Site, y = Normalized.Abundance.rplB, fill = Method)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Gene) +
  theme_classic() +
  ylab("rplB-normalized abundance") +
  theme(axis.text.x = element_text(angle = 90, size = 10, 
                                   hjust=0.95,vjust=0.2))

ggplot(mixed, aes(x = Method, y = Normalized.Abundance.rplB)) +
  geom_boxplot() +
  geom_jitter(aes(color = Site), width = 0.3, size = 3) +
  facet_wrap(~Gene) +
  theme_classic() +
  scale_color_brewer(palette = "Paired") +
  ylab("rplB-normalized abundance") +
  theme(axis.text.x = element_text(angle = 90, size = 10, 
                                   hjust=0.95,vjust=0.2))

####################
#look at rplB comps#
####################

#summarise xander rplB
rplB.xander.summary <- rplB.xander %>%
  group_by(Site) %>%
  summarise(rplB = sum(Abundance)) %>%
  mutate(Method = "Xander")

#add method to cog rplb
cogs.rplB.summary <- rplB %>%
  mutate(Method = "COGs") %>%
  select(Site, rplB, Method)

#change site from c to cen
cogs.rplB.summary$Site <- gsub("C", "Cen", cogs.rplB.summary$Site)

#bind datasets
rplB.summary <- rbind(cogs.rplB.summary, rplB.xander.summary)

#plot data
ggplot(rplB.summary, aes(x = Site, y = rplB, fill = Method)) +
  geom_col(position = "dodge", color = "black") +
  theme_classic() +
  ylab("Gene count") +
  theme(axis.text.x = element_text(angle = 90, size = 10, 
                                   hjust=0.95,vjust=0.2))
