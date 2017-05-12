#load dependencies 
library(phyloseq)
library(ape)
library(biomformat)
library(vegan)
library(tidyverse)
library(reshape2)

#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis
wd=print(getwd())

##################################
#PHYLUM_LEVEL_RESPONSES_WITH_RPLB#
##################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names=list.files(pattern="*rplB_45_taxonabund.txt")
data <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read.csv(X))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#remove rows that include lineage match name (only need taxon)
#aka remove duplicate data
data <- data[!grepl(";", data$Taxon.Abundance.Fraction.Abundance),]
data <- data[!grepl("Lineage", data$Taxon.Abundance.Fraction.Abundance),]

#split columns 
data <- data %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = Taxon.Abundance.Fraction.Abundance, 
           into = c("Taxon", "Abundance", "Fraction.Abundance"), 
           sep = "\t") %>%
  select(-junk) %>%
  group_by(Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data$Fraction.Abundance <- as.numeric(data$Fraction.Abundance)
data$Abundance <- as.numeric(data$Abundance)

#change site from "cen" to "Cen" so it matches metadata
data$Site <- gsub("cen", "Cen", data$Site)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised <- data %>%
  summarise(Total = sum(Fraction.Abundance), rplB = sum(Abundance))

#change site from "cen" to "Cen" so it matches metadata
summarised$Site <- gsub("cen", "Cen", summarised$Site)

#decast for abundance check
dcast=acast(data, Taxon ~ Site, value.var = "Fraction.Abundance")

#call na's zeros
dcast[is.na(dcast)] =0

#order based on abundance
order.dcast=dcast[order(rowSums(dcast),decreasing=TRUE),]

#melt data
melt=melt(order.dcast,id.vars=row.names(order.dcast), variable.name= "Site", value.name="Fraction.Abundance" )

#adjust colnames of melt
colnames(melt)=c("Taxon", "Site", "Fraction.Abundance")

#read in metadata
meta=read.delim(file = paste(wd, "/data/Centralia_full_map.txt", sep=""), sep=" ")

#join metadata with regular data
joined=inner_join(melt, meta, by="Site")

#average by fire history
#group data
grouped=group_by(joined, Taxon, Classification)

#calculate
history=summarise(grouped, N=length(Fraction.Abundance), Average=mean(Fraction.Abundance))

#plot
(phylum.plot=(ggplot(history, aes(x=Taxon, y=Average)) +
                geom_point(size=2) +
                facet_wrap(~Classification, ncol = 1) +
                labs(x="Phylum", y="Mean relative abundance")+
                theme(axis.text.x = element_text(angle = 90, size = 12, hjust=0.95,vjust=0.2))))

#save plot
ggsave(phylum.plot, filename = paste(wd, "/figures/phylum.responses.png", sep=""), width = 5, height = 5)

#########################
#RplB DIVERSITY ANALYSIS#
#########################
#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#remove Cen16 from metadata since we don't have rplB info yet
meta=meta[!grepl("Cen16", meta$Site),]

#read in distance matrix
rplB=read.delim(file = paste(wd, "/data/rformat_dist_0.03.txt", sep=""))

#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#add row names back
rownames(rplB)=rplB[,1]

#remove first column
rplB=rplB[,-1]

#make data matrix
rplB=data.matrix(rplB)

#remove first column
rplB=rplB[,-1]

#make an output of total gene count per site
rplB.gcounts=rowSums(rplB)

#otu table
otu=otu_table(rplB, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = FALSE)

#rarefy
rare=rarefy_even_depth(otu, sample.size = min(sample_sums(otu)), rngseed = TRUE)

#check curve
rarecurve(rare, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo=merge_phyloseq(rare, metad)

#plot phylo richness
(richness=plot_richness(phylo, x="Sample", shape="Classification", color = "SoilTemperature_to10cm"))

#save plot 
ggsave(richness, filename = paste(wd, "/figures/richness.png", sep=""), 
       width = 15)

#calculate evenness
s=specnumber(rare)
h=diversity(rare, index="shannon")
plieou=h/log(s)

#save evenness number
write.table(plieou, file = paste(wd, "/output/evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou=data.frame(plieou)

#add site column to evenness data
plieou$Site=rownames(plieou)

#merge evenness information with fire classification
plieou=inner_join(plieou, meta)

#make color pallette
GnYlOrRd=colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

#plot evenness by fire classification
(evenness <- ggplot(plieou, aes(x = Classification, y = plieou)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3, width = 0.2) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(evenness, filename = paste(wd, "/figures/evenness.png", sep=""), 
       width = 7, height = 5)

#plot Bray Curtis ordination
ord <- ordinate(phylo, method="PCoA", distance="bray")
(bc.ord=plot_ordination(phylo, ord, shape="Classification", title="Bray Curtis") +
  geom_point(aes(color = SoilTemperature_to10cm), size=5) +
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
  theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord, filename = paste(wd, "/figures/braycurtis.ord.png", sep=""), 
       width = 6, height = 5)


#plot Sorenson ordination
ord.sor <- ordinate(phylo, method="PCoA", distance="bray", binary = TRUE)
(sorenson.ord=plot_ordination(phylo, ord.sor, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord, filename = paste(wd, "/figures/sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#make object phylo with tree and biom info
tree <- read.tree(file = paste(wd, "/data/rplB_0.03_tree.nwk", sep=""))
tree <- phy_tree(tree)

#merge
phylo=merge_phyloseq(tree, rare, metad)

#plot unweighted Unifrac ordination
uni.u.ord.rplB <- ordinate(phylo, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.rplB=plot_ordination(phylo, uni.u.ord.rplB, shape="Classification", 
                                title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.rplB, filename = paste(wd, "/figures/rplB.u.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.rplB <- ordinate(phylo, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.rplB=plot_ordination(phylo, uni.w.ord.rplB, shape="Classification", 
                                title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.rplB, filename = paste(wd, "/figures/rplB.w.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot tree
(tree.plot <- plot_tree(phylo, shape = "Classification", size = "abundance",
                        color = "SoilTemperature_to10cm", label.tips=NULL, 
                        text.size=2, ladderize="left", base.spacing = 0.04) +
    theme(legend.position = "right", legend.title = element_text(size=16),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(tree.plot, filename = paste(wd, "/figures/rplb.tree.png", sep=""), 
       width = 10, height = 40)

#########################
#ArsB DIVERSITY ANALYSIS#
#########################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#remove Cen16 from metadata since we don't have arsB info yet
meta=meta[which(meta$Site == "Cen17" | meta$Site == "Cen10"),]
meta$Site <- as.character(meta$Site)

#read in distance matrix
arsB=read.delim(file = paste(wd, "/data/arsB_rformat_dist_0.03.txt", sep=""))

#add row names back
rownames(arsB)=arsB[,1]

#remove first column
arsB=arsB[,-1]

#make data matrix
arsB=data.matrix(arsB)

#remove first column
arsB=arsB[,-1]

#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#otu table
otu=otu_table(arsB, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = FALSE)

#rarefy
rare=rarefy_even_depth(otu, sample.size = min(sample_sums(otu)), rngseed = TRUE)

#check curve
rarecurve(rare, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo=merge_phyloseq(rare, metad)

#plot phylo richness
(richness=plot_richness(phylo, x="Sample", shape="Classification", color = "SoilTemperature_to10cm"))

#save plot 
ggsave(richness, filename = paste(wd, "/figures/arsB.richness.png", sep=""), 
       width = 15)

#calculate evenness
s=specnumber(rare)
h=diversity(rare, index="shannon")
plieou=h/log(s)

#save evenness number
write.table(plieou, file = paste(wd, "/output/arsB.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou=data.frame(plieou)

#add site column to evenness data
plieou$Site=rownames(plieou)

#merge evenness information with fire classification
plieou=inner_join(plieou, meta)

#make color pallette
GnYlOrRd=colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

#plot evenness by fire classification
(evenness <- ggplot(plieou, aes(x = Classification, y = plieou)) +
    geom_point(aes(color = SoilTemperature_to10cm), size=3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(evenness, filename = paste(wd, "/figures/arsB.evenness.png", sep=""), 
       width = 7, height = 5)

#make an output of total gene count per site
arsB.gcounts=rowSums(arsB)

#make object phylo with tree and biom info
tree <- read.tree(file = paste(wd, "/data/arsB_0.03_tree.nwk", sep=""))
tree <- phy_tree(tree)

#merge
phylo=merge_phyloseq(tree, rare, metad)

#plot tree
(arsB.tree.plot <- plot_tree(phylo, color = "SoilTemperature_to10cm", size = "abundance",
                        shape = "Classification", label.tips=NULL, 
                        text.size=2, ladderize="left", base.spacing = 0.03, sizebase = 2) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(arsB.tree.plot, filename = paste(wd, "/figures/arsB.tree.png", sep=""))

##############################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR ARSB#
##############################################
#read in each file (one per site)
cen10 <- read_delim(file = paste(wd, "/data/arsB_taxonabund_cen10.txt", sep = ""), 
                    col_names = TRUE, delim = "\t")
cen17 <- read_delim(file = paste(wd, "/data/arsB_taxonabund_cen17.txt", sep = ""), 
                    col_names = TRUE, delim = "\t")

#make column for organism name
cen10 <- cen10 %>%
  mutate(Site = "cen10", Gene = "arsB", Census = 8744.20, rplB = 3040) %>%
  mutate(Normalized.Abundance.rplB = Abundance / rplB, 
         Normalized.Abundance.census = Abundance / Census) %>%
  separate(Taxon, into = c("coded_by", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), sep = ",definition=")

cen17 <- cen17 %>%
  mutate(Site = "cen17", Gene = "arsB", Census = 3867.45, rplB = 4928) %>%
  mutate(Normalized.Abundance.rplB = Abundance / rplB, 
         Normalized.Abundance.census = Abundance / Census) %>%  
  separate(Taxon, into = c("coded_by", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), sep = ",definition=")

#join data together
arsB <- rbind(cen17, cen10)
arsB.high <- arsB[which(arsB$Abundance > 10),]

#plot data
(arsB.abundance.plot <- ggplot(arsB, aes(x = Site, y = Normalized.Abundance.rplB)) +
  geom_bar(stat = "identity") +
    ylab("arsB Abundance (normalized to rplB)"))

ggsave(arsB.abundance.plot, filename = paste(wd, "/figures/arsB.abundance.png", sep=""))

ggplot(arsB, aes(x = Site, y = Normalized.Abundance.census)) +
  geom_bar(stat = "identity")

ggplot(arsB, aes(x = organism, y = Normalized.Abundance.census)) +
  geom_point(aes(color = Site)) +
  coord_flip()

#plot 
(arsB.taxon.plot <- ggplot(arsB, aes(x = organism, y = Normalized.Abundance.rplB)) +
  geom_point(aes(color = Site)) +
  ylab("arsB Abundance (normalized to rplB)") +
  xlab("Taxon") +
  coord_flip())

ggsave(arsB.taxon.plot, filename = paste(wd, "/figures/arsB.abundance.taxon.png", sep=""), 
       width = 7, height = 5)

#########################
#ACR3 DIVERSITY ANALYSIS#
#########################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#remove Cen16 from metadata since we don't have acr3 info yet
meta=meta[which(meta$Site == "Cen14" | meta$Site == "Cen03" | meta$Site == "Cen01" | meta$Site == "Cen04" | meta$Site == "Cen07" | meta$Site == "Cen15" | meta$Site == "Cen12"),]
meta$Site <- as.character(meta$Site)

#read in distance matrix
acr3=read.delim(file = paste(wd, "/data/acr3_rformat_dist_0.03.txt", sep=""))

#add row names back
rownames(acr3)=acr3[,1]

#remove first column
acr3=acr3[,-1]

#make data matrix
acr3=data.matrix(acr3)

#remove first column
acr3=acr3[,-1]

#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#otu table
otu=otu_table(acr3, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = FALSE)

#rarefy
rare=rarefy_even_depth(otu, sample.size = min(sample_sums(otu)), rngseed = TRUE)

#check curve
rarecurve(rare, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo=merge_phyloseq(rare, metad)

#plot phylo richness
(richness=plot_richness(phylo, x="Sample", shape="Classification", color = "SoilTemperature_to10cm"))

#save plot 
ggsave(richness, filename = paste(wd, "/figures/acr3.richness.png", sep=""), 
       width = 15)

#calculate evenness
s=specnumber(rare)
h=diversity(rare, index="shannon")
plieou=h/log(s)

#save evenness number
write.table(plieou, file = paste(wd, "/output/acr3.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou=data.frame(plieou)

#add site column to evenness data
plieou$Site=rownames(plieou)

#merge evenness information with fire classification
plieou=inner_join(plieou, meta)

#make color pallette
GnYlOrRd=colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

#plot evenness by fire classification
(evenness <- ggplot(plieou, aes(x = Classification, y = plieou)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(evenness, filename = paste(wd, "/figures/acr3.evenness.png", sep=""), 
       width = 7, height = 5)

#make an output of total gene count per site
acr3.gcounts=rowSums(acr3)

#make object phylo with tree and biom info
tree <- read.tree(file = paste(wd, "/data/acr3_0.03_tree.nwk", sep=""))
tree <- phy_tree(tree)

#merge
phylo=merge_phyloseq(tree, rare, metad)


#plot Bray Curtis ordination
ord.acr3 <- ordinate(phylo, method="PCoA", distance="bray")
(bc.ord.acr3=plot_ordination(phylo, ord.acr3, shape="Classification", title="Bray Curtis") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord.acr3, filename = paste(wd, "/figures/acr3.braycurtis.ord.png", sep=""), 
       width = 6, height = 5)

#plot Sorenson ordination
s.ord.acr3 <- ordinate(phylo, method="PCoA", distance="bray", binary=TRUE)
(sorenson.ord.acr3=plot_ordination(phylo, s.ord.acr3, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord.acr3, filename = paste(wd, "/figures/acr3.sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#plot unweighted Unifrac ordination
uni.u.ord.acr3 <- ordinate(phylo, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.acr3=plot_ordination(phylo, uni.u.ord.acr3, shape="Classification", 
                              title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.acr3, filename = paste(wd, "/figures/acr3.u.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.acr3 <- ordinate(phylo, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.acr3=plot_ordination(phylo, uni.w.ord.acr3, shape="Classification", 
                              title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.acr3, filename = paste(wd, "/figures/acr3.w.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot tree
(acr3.tree.plot <- plot_tree(phylo, color = "SoilTemperature_to10cm", size = "abundance",
                        shape = "Classification", label.tips=NULL, 
                        text.size=2, ladderize="left", base.spacing = 0.04) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(acr3.tree.plot, filename = paste(wd, "/figures/acr3.tree.png", sep=""), 
       height = 10, width = 10)


##############################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR ACR3#
##############################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names.acr3=list.files(pattern="*acr3_45_taxonabund.txt")
data.acr3 <- do.call(rbind, lapply(names.acr3, function(X) {
  data.frame(id = basename(X), read_table(X))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#remove NA rows 
data.acr3 <- data.acr3[!is.na(data.acr3$Taxon.Abundance.Fraction.Abundance),]

#split columns 
data.acr3 <- data.acr3 %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = Taxon.Abundance.Fraction.Abundance, 
           into = c("Taxon", "Abundance", "Fraction.Abundance"), 
           sep = "\t") %>%
  select(-junk) %>%
  group_by(Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.acr3$Fraction.Abundance <- as.numeric(data.acr3$Fraction.Abundance)
data.acr3$Abundance <- as.numeric(data.acr3$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.acr3 <- data.acr3 %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#change "cen" to "Cen" so it matches outside data
data.acr3$Site <- gsub("cen", "Cen", data.acr3$Site)

#read in microbe census data
census <- read_delim(file = paste(wd, "/data/microbe_census.txt", sep = ""), delim = "\t", 
                     col_types = list(col_character(), col_number(), col_number(), col_number()))

#read in arsenic and temperature data
meta <- data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#make column for organism name and join with microbe census data and normalize to it
data.acr3 <- data.acr3 %>%
  left_join(census, by = "Site") %>%
  left_join(summarised, by = "Site") %>%
  left_join(meta, by = "Site") %>%
  separate(Taxon, into = c("coded_by", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), sep = ",definition=") %>%
  rename(Temp = SoilTemperature_to10cm) %>%
  select(Site, As_ppm, Temp, Classification, organism, Abundance, 
         Fraction.Abundance, GE, rplB) %>%
  mutate(Normalized.Abundance.census = Abundance / GE, 
         Normalized.Abundance.rplB = Abundance / rplB)

###TEMPORARY
#remove Cen01 since it cant be normalized to rplB yet
data.acr3 <- data.acr3[-which(data.acr3$Site == "Cen01"),]

#plot data
(acr3.abundance.plot <- ggplot(data.acr3, aes(x = Site, y = Normalized.Abundance.rplB, 
                                              fill = Classification)) +
  geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow")) +
    ylab("acr3 Abundance (normalized to rplB)") +
    theme_classic())

ggsave(acr3.abundance.plot, filename = paste(wd, "/figures/acr3.abundance.png", sep=""))

(acr3.abundance.census.plot <- ggplot(data.acr3, aes(x = Site, 
                                                     y = Normalized.Abundance.census, 
                                                     fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow")) +
    ylab("acr3 Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(acr3.abundance.census.plot, filename = paste(wd, "/figures/acr3.abundance.census.png",
                                                    sep=""))


(acr3.abundance.census.taxon.plot <- ggplot(data.acr3, aes(x = organism, 
                                                    y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("acr3 abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(acr3.abundance.census.taxon.plot, 
       filename = paste(wd, "/figures/acr3.abundance.census.taxon.png", sep=""), height = 6.5)

(acr3.abundance.taxon.plot <- ggplot(data.acr3, aes(x = organism, 
                                               y = Normalized.Abundance.rplB)) +
  geom_point(aes(color = Temp, shape = Classification)) +
    ylab("acr3 abundance (normalized to rplB)") +
    xlab("Taxon") +
  coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(acr3.abundance.taxon.plot, 
       filename = paste(wd, "/figures/acr3.abundance.taxon.png", sep=""), height = 6.5)

#########################
#AIOA DIVERSITY ANALYSIS#
#########################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#remove Cen16 from metadata since we don't have aioA info yet
meta=meta[which(meta$Site == "Cen07" | meta$Site == "Cen10"),]
meta$Site <- as.character(meta$Site)

#read in distance matrix
aioA=read.delim(file = paste(wd, "/data/aioA_rformat_dist_0.03.txt", sep=""))

#add row names back
rownames(aioA)=aioA[,1]

#remove first column
aioA=aioA[,-1]

#make data matrix
aioA=data.matrix(aioA)

#remove first column
aioA=aioA[,-1]


#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#otu table
otu=otu_table(aioA, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = FALSE)

#rarefy
rare=rarefy_even_depth(otu, sample.size = min(sample_sums(otu)), rngseed = TRUE)

#check curve
rarecurve(rare, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo=merge_phyloseq(rare, metad)

#plot phylo richness
(richness=plot_richness(phylo, x="Sample", shape="Classification", color = "SoilTemperature_to10cm"))

#save plot 
ggsave(richness, filename = paste(wd, "/figures/aioA.richness.png", sep=""), 
       width = 15)

#calculate evenness
s=specnumber(rare)
h=diversity(rare, index="shannon")
plieou=h/log(s)

#save evenness number
write.table(plieou, file = paste(wd, "/output/aioA.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou=data.frame(plieou)

#add site column to evenness data
plieou$Site=rownames(plieou)

#merge evenness information with fire classification
plieou=inner_join(plieou, meta)

#plot evenness by fire classification
(evenness <- ggplot(plieou, aes(x = Classification, y = plieou)) +
    geom_point(aes(color = SoilTemperature_to10cm), size=3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(evenness, filename = paste(wd, "/figures/aioA.evenness.png", sep=""), 
       width = 7, height = 5)

#make an output of total gene count per site
aioA.gcounts=rowSums(aioA)

#make object phylo with tree and biom info
tree <- read.tree(file = paste(wd, "/data/aioA_0.03_tree.nwk", sep=""))
tree <- phy_tree(tree)

#merge
phylo=merge_phyloseq(tree, rare, metad)

#plot tree
(aioA.tree.plot <- plot_tree(phylo, color = "SoilTemperature_to10cm", size = "abundance",
                        shape = "Classification", label.tips="taxa_names", 
                        text.size=2, ladderize="left", base.spacing = 0.03, sizebase = 2) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(aioA.tree.plot, filename = paste(wd, "/figures/aioA.tree.png", sep=""))


##############################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR AIOA#
##############################################

#read in each file (one per site)
cen10 <- read_delim(file = paste(wd, "/data/aioA_taxonabund_cen10.txt", sep = ""), 
                    col_names = TRUE, delim = "\t")
cen07 <- read_delim(file = paste(wd, "/data/aioA_taxonabund_cen07.txt", sep = ""), 
                    col_names = TRUE, delim = "\t")

#make column for organism name
cen10 <- cen10 %>%
  mutate(Site = "cen10", Gene = "aioA", Census = 8744.20, rplB = 3040) %>%
  mutate(Normalized.Abundance.rplB = Abundance / rplB, 
         Normalized.Abundance.census = Abundance / Census) %>%
  separate(Taxon, into = c("coded_by", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), sep = ",definition=")

cen07 <- cen07 %>%
  mutate(Site = "cen07", Gene = "aioA", Census = 3584.52, rplB = 1347.908) %>%
  mutate(Normalized.Abundance.rplB = Abundance / rplB, 
         Normalized.Abundance.census = Abundance / Census) %>%  
  separate(Taxon, into = c("coded_by", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), sep = ",definition=")

#join data together
aioA <- rbind(cen07, cen10)
aioA.high <- aioA[which(aioA$Abundance > 10),]

#plot data
(aioA.abundance.plot <- ggplot(aioA, aes(x = Site, y = Normalized.Abundance.rplB)) +
  geom_bar(stat = "identity") +
    ylab("aioA Abundance (normalized to rplB)"))

ggsave(aioA.abundance.plot, filename = paste(wd, "/figures/aioA.abundance.png", sep = ""))

ggplot(aioA, aes(x = Site, y = Normalized.Abundance.census)) +
  geom_bar(stat = "identity")

ggplot(aioA, aes(x = organism, y = Normalized.Abundance.census)) +
  geom_point(aes(color = Site)) +
  coord_flip()

(aioA.abundance.taxon.plot <- ggplot(aioA, aes(x = organism, 
                                               y = Normalized.Abundance.rplB)) +
  geom_point(aes(color = Site)) +
    ylab("aioA abundance (normalized to rplB)") +
    xlab("Taxon") +
  coord_flip())
ggsave(aioA.abundance.taxon.plot, filename = paste(wd, "/figures/aioA.abundance.taxon.png", sep = ""))















