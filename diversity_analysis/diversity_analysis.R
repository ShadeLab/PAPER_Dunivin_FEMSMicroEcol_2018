#load dependencies 
library(phyloseq)
library(ape)
library(biomformat)
library(vegan)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(taxize)


#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis
wd <- print(getwd())

#make color pallette
GnYlOrRd <- colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

#read in microbe census data
census <- read_delim(file = paste(wd, "/data/microbe_census.txt", sep = ""),
                     delim = "\t", col_types = list(col_character(), col_number(),
                                                    col_number(), col_number()))

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
(richness=plot_richness(phylo, x="Classification", color = "SoilTemperature_to10cm"))

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
otu.arsB=otu_table(arsB, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu.arsB, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = FALSE)

#rarefy
rare.arsB=rarefy_even_depth(otu.arsB, sample.size = min(sample_sums(otu.arsB)), rngseed = TRUE)

#check curve
rarecurve(rare.arsB, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.arsB=merge_phyloseq(rare.arsB, metad)

#plot phylo richness
(richness=plot_richness(phylo.arsB, x="Classification", color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(richness, filename = paste(wd, "/figures/arsB.richness.png", sep=""), 
       width = 15)

#calculate evenness
s.arsB=specnumber(rare.arsB)
h.arsB=diversity(rare.arsB, index="shannon")
plieou.arsB=h.arsB/log(s.arsB)

#save evenness number
write.table(plieou.arsB, file = paste(wd, "/output/arsB.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou.arsB=data.frame(plieou.arsB)

#add site column to evenness data
plieou.arsB$Site=rownames(plieou.arsB)

#merge evenness information with fire classification
plieou.arsB=inner_join(plieou.arsB, meta)

#plot evenness by fire classification
(evenness <- ggplot(plieou.arsB, aes(x = Classification, y = plieou.arsB)) +
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
tree.arsB <- read.tree(file = paste(wd, "/data/arsB_0.03_tree.nwk", sep=""))
tree.arsB <- phy_tree(tree.arsB)

#merge
phylo.arsB=merge_phyloseq(tree.arsB, rare.arsB, metad)

#plot tree
(arsB.tree.plot <- plot_tree(phylo.arsB, color = "SoilTemperature_to10cm", size = "abundance",
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
meta=meta[-which(meta$Site == "Cen10"),]
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
otu.acr3=otu_table(acr3, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu.acr3, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = FALSE)

#rarefy
rare.acr3=rarefy_even_depth(otu.acr3, sample.size = min(sample_sums(otu.acr3)), rngseed = TRUE)

#check curve
rarecurve(rare.acr3, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.acr3=merge_phyloseq(rare.acr3, metad)

#plot phylo richness
(richness=plot_richness(phylo.acr3, x="Classification", color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(richness, filename = paste(wd, "/figures/acr3.richness.png", sep=""), 
       width = 15)

#calculate evenness
s.acr3=specnumber(rare.acr3)
h.acr3=diversity(rare.acr3, index="shannon")
plieou.acr3=h.acr3/log(s.acr3)

#save evenness number
write.table(plieou.acr3, file = paste(wd, "/output/acr3.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou.acr3=data.frame(plieou.acr3)

#add site column to evenness data
plieou.acr3$Site=rownames(plieou.acr3)

#merge evenness information with fire classification
plieou.acr3=inner_join(plieou.acr3, meta)

#plot evenness by fire classification
(evenness <- ggplot(plieou.acr3, aes(x = Classification, y = plieou.acr3)) +
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
tree.acr3 <- read.tree(file = paste(wd, "/data/acr3_0.03_tree.nwk", sep=""))
tree.acr3 <- phy_tree(tree.acr3)

#merge
phylo.acr3=merge_phyloseq(tree.acr3, rare.acr3, metad)


#plot Bray Curtis ordination
ord.acr3 <- ordinate(phylo.acr3, method="PCoA", distance="bray")
(bc.ord.acr3=plot_ordination(phylo.acr3, ord.acr3, shape="Classification", title="Bray Curtis") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord.acr3, filename = paste(wd, "/figures/acr3.braycurtis.ord.png", sep=""), 
       width = 6, height = 5)

#plot Sorenson ordination
s.ord.acr3 <- ordinate(phylo.acr3, method="PCoA", distance="bray", binary=TRUE)
(sorenson.ord.acr3=plot_ordination(phylo.acr3, s.ord.acr3, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord.acr3, filename = paste(wd, "/figures/acr3.sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#plot unweighted Unifrac ordination
uni.u.ord.acr3 <- ordinate(phylo.acr3, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.acr3=plot_ordination(phylo.acr3, uni.u.ord.acr3, shape="Classification", 
                              title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.acr3, filename = paste(wd, "/figures/acr3.u.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.acr3 <- ordinate(phylo.acr3, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.acr3=plot_ordination(phylo.acr3, uni.w.ord.acr3, shape="Classification", 
                              title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.acr3, filename = paste(wd, "/figures/acr3.w.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot tree
(acr3.tree.plot <- plot_tree(phylo.acr3, color = "SoilTemperature_to10cm", 
                             size = "abundance",
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

#get taxonomy for organisms
ncbi <- tax_name(query = data.acr3$organism, get = c("genus", "class", "phylum"), db = "ncbi")

#change query column to "organism"
ncbi$organism <- ncbi$query

#join taxanomic information with abundance information
data.acr3.t <- data.acr3 %>%
  left_join(ncbi, by = "organism") %>%
  unique()

#change NA phylum to metagenome
data.acr3.t$phylum[is.na(data.acr3.t$phylum)] = "Metagenome"

#order based on temperature
data.acr3.t$Site <- factor(data.acr3.t$Site, 
                           levels = data.acr3.t$Site[order(data.acr3.t$Temp)])

#prep colors
n <- 32
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color.acr3 <- print(sample(col_vector, n))


(acr3.abundance.tax.bar.plot <- ggplot(data.acr3.t, aes(x = Site, y = Normalized.Abundance.census, fill = phylum)) +
    geom_bar(stat = "identity") +
    ylab("acr3 Abundance (normalized to genome equivalents)") +
    scale_fill_manual(values = color.acr3) +
    theme_classic(base_size = 12))
ggsave(acr3.abundance.tax.bar.plot, filename = paste(wd, "/figures/acr3.phylum.eps", sep=""), width = 10)

#order based on temperature
data.acr3$Site <- factor(data.acr3$Site, 
                           levels = data.acr3$Site[order(data.acr3$Temp)])

#plot data
(acr3.abundance.plot <- ggplot(data.acr3, aes(x = Site, y = Normalized.Abundance.rplB, 
                                              fill = Classification)) +
  geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow", "green")) +
    ylab("acr3 Abundance (normalized to rplB)") +
    theme_classic())

ggsave(acr3.abundance.plot, filename = paste(wd, "/figures/acr3.abundance.png", sep=""))

(acr3.abundance.census.plot <- ggplot(data.acr3, aes(x = Site, 
                                                     y = Normalized.Abundance.census, 
                                                     fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("firebrick2", "yellow1", "green3")) +
    ylab("acr3 Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(acr3.abundance.census.plot, 
       filename = paste(wd, "/figures/acr3.abundance.census.png", sep=""))


(acr3.abundance.census.taxon.plot <- ggplot(data.acr3, aes(x = organism, 
                                                    y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("acr3 abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(acr3.abundance.census.taxon.plot, 
       filename = paste(wd, "/figures/acr3.abundance.census.taxon.png", sep=""), height = 10)

(acr3.abundance.taxon.plot <- ggplot(data.acr3, aes(x = organism, 
                                               y = Normalized.Abundance.rplB)) +
  geom_point(aes(color = Temp, shape = Classification)) +
    ylab("acr3 abundance (normalized to rplB)") +
    xlab("Taxon") +
  coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(acr3.abundance.taxon.plot, 
       filename = paste(wd, "/figures/acr3.abundance.taxon.png", sep=""), height = 10)

#summarise acr3 information
data.acr3.sum <- data.acr3 %>%
  group_by(Classification, Site) %>%
  summarise(Abundance.census.tot = sum(Normalized.Abundance.census), 
            Abundance.rplB.tot = sum(Normalized.Abundance.rplB))

(acr3.abundance.census.box <- ggplot(data.acr3.sum, aes(x = Classification, 
                                                    y = Abundance.census.tot)) +
    geom_boxplot() +
    geom_jitter(width = 0.3) +
    ylab("acr3 abundance (normalized to genome equivalents)"))

ggsave(acr3.abundance.census.box, 
       filename = paste(wd, "/figures/acr3.abundance.census.box.png", sep=""))

#########################
#AIOA DIVERSITY ANALYSIS#
#########################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#remove Cen16 from metadata since we don't have aioA info yet
meta=meta[-which(meta$Site == "Cen14" | meta$Site == "Cen15" | meta$Site == "Cen16"),]
meta$Site <- as.character(meta$Site)

#read in distance matrix
aioA=read.delim(file = paste(wd, "/data/aioA_rformat_dist_0.03.txt", sep=""))

##temporary: remove Cen12 since it only has 1 hit
aioA <- aioA[-which(aioA$X == "Cen12" | aioA$X == "Cen17"),]

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
otu.aioA=otu_table(aioA, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu.aioA, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

#rarefy (not for now since we do not have enough data?)
rare.aioA <- otu.aioA
#rare.aioA=rarefy_even_depth(otu.aioA, sample.size = median(sample_sums(otu.aioA)), rngseed = TRUE)

#check curve
rarecurve(rare.aioA, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.aioA=merge_phyloseq(rare.aioA, metad)

#plot phylo richness
(richness=plot_richness(phylo.aioA, x="Classification", 
                        color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(richness, filename = paste(wd, "/figures/aioA.richness.png", sep=""), 
       width = 15)

#calculate evenness
s.aioA=specnumber(rare.aioA)
h.aioA=diversity(rare.aioA, index="shannon")
plieou.aioA=h.aioA/log(s.aioA)

#save evenness number
write.table(plieou.aioA, file = paste(wd, "/output/aioA.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou.aioA=data.frame(plieou.aioA)

#add site column to evenness data
plieou.aioA$Site=rownames(plieou.aioA)

#merge evenness information with fire classification
plieou.aioA=inner_join(plieou.aioA, meta)

#plot evenness by fire classification
(evenness <- ggplot(plieou.aioA, aes(x = Classification, y = plieou.aioA)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3, width = 0.15) +
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
tree.aioA <- read.tree(file = paste(wd, "/data/aioA_0.03_tree.nwk", sep=""))
tree.aioA <- phy_tree(tree.aioA)

#merge
phylo.aioA=merge_phyloseq(tree.aioA, rare.aioA, metad)

#plot tree
(aioA.tree.plot <- plot_tree(phylo.aioA, color = "SoilTemperature_to10cm", 
                             size = "abundance",
                        shape = "Classification", label.tips="taxa_names", 
                        text.size=2, ladderize="left", base.spacing = 0.03, sizebase = 2) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(aioA.tree.plot, filename = paste(wd, "/figures/aioA.tree.png", sep=""))

#plot Bray Curtis ordination
ord.aioA <- ordinate(phylo.aioA, method="PCoA", distance="bray")
(bc.ord.aioA=plot_ordination(phylo.aioA, ord.aioA, shape="Classification",
                             title="Bray Curtis") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord.aioA, filename = paste(wd, "/figures/aioA.braycurtis.ord.png", sep=""), 
       width = 6, height = 5)

#plot Sorenson ordination
s.ord.aioA <- ordinate(phylo.aioA, method="PCoA", distance="bray", binary=TRUE)
(sorenson.ord.aioA=plot_ordination(phylo.aioA, s.ord.aioA, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord.aioA, filename = paste(wd, "/figures/aioA.sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#plot unweighted Unifrac ordination
uni.u.ord.aioA <- ordinate(phylo.aioA, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.aioA=plot_ordination(phylo.aioA, uni.u.ord.aioA, shape="Classification", 
                                title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.aioA, filename = paste(wd, "/figures/aioA.u.unifrac.ord.png", 
                                        sep=""), width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.aioA <- ordinate(phylo.aioA, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.aioA=plot_ordination(phylo.aioA, uni.w.ord.aioA, 
                                shape="Classification", title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.aioA, filename = paste(wd, "/figures/aioA.w.unifrac.ord.png", 
                                        sep=""), width = 6, height = 5)


##############################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR AIOA#
##############################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names.aioA=list.files(pattern="*aioA_45_taxonabund.txt")
data.aioA <- do.call(rbind, lapply(names.aioA, function(X) {
  data.frame(id = basename(X), read_table(X))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#split columns 
data.aioA <- data.aioA %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = Taxon.Abundance.FractionAbundance, 
           into = c("Taxon", "Abundance", "Fraction.Abundance"), 
           sep = "\t") %>%
  select(-junk) %>%
  group_by(Site)

#manually fill in cen 12
data.aioA <- data.aioA[-which(data.aioA$Site == "cen12"),]

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.aioA$Fraction.Abundance <- as.numeric(data.aioA$Fraction.Abundance)
data.aioA$Abundance <- as.numeric(data.aioA$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.aioA <- data.aioA %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#change "cen" to "Cen" so it matches outside data
data.aioA$Site <- gsub("cen", "Cen", data.aioA$Site)

#read in arsenic and temperature data
meta <- data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#make column for organism name and join with microbe census data and normalize to it
data.aioA <- data.aioA %>%
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

#get taxonomy for organisms
ncbi <- tax_name(query = data.aioA$organism, get = c("genus", "class", "phylum"), db = "ncbi")

#change query column to "organism"
ncbi$organism <- ncbi$query

#join taxanomic information with abundance information
data.aioA.t <- data.aioA %>%
  left_join(ncbi, by = "organism") %>%
  unique()

#prep colors
n <- 19
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color <- print(sample(col_vector, n))

#order based on temperature
data.aioA.t$Site <- factor(data.aioA.t$Site, 
                           levels = data.aioA.t$Site[order(data.aioA.t$Temp)])
(aioA.abundance.tax.bar.plot <- ggplot(data.aioA.t, aes(x = Site, y = Normalized.Abundance.census, fill = genus)) +
    geom_bar(stat = "identity") +
    ylab("acr3 Abundance (normalized to genome equivalents)") +
    scale_fill_manual(values = color) +
    theme_classic(base_size = 12))
ggsave(aioA.abundance.tax.bar.plot, filename = paste(wd, "/figures/aioA.genus.eps", sep=""), width = 10)


#plot data
(aioA.abundance.plot <- ggplot(data.aioA, aes(x = Site, y = Normalized.Abundance.rplB, 
                                              fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow", "green")) +
    ylab("aioA Abundance (normalized to rplB)") +
    theme_classic())

ggsave(aioA.abundance.plot, filename = paste(wd, "/figures/aioA.abundance.png", sep=""))

(aioA.abundance.census.plot <- ggplot(data.aioA, aes(x = Site, 
                                                     y = Normalized.Abundance.census, 
                                                     fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow", "green")) +
    ylab("aioA Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(aioA.abundance.census.plot, filename = paste(wd, "/figures/aioA.abundance.census.png",
                                                    sep=""))


(aioA.abundance.census.taxon.plot <- ggplot(data.aioA, aes(x = organism, 
                                                           y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("aioA abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(aioA.abundance.census.taxon.plot, 
       filename = paste(wd, "/figures/aioA.abundance.census.taxon.png", sep=""), height = 10)

(aioA.abundance.taxon.plot <- ggplot(data.aioA, aes(x = organism, 
                                                    y = Normalized.Abundance.rplB)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("aioA abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(aioA.abundance.taxon.plot, 
       filename = paste(wd, "/figures/aioA.abundance.taxon.png", sep=""), height = 10)

#########################
#arsM DIVERSITY ANALYSIS#
#########################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#remove Cen16 from metadata since we don't have arsM info yet
meta=meta[-which(meta$Site == "Cen05" | meta$Site == "Cen10" | meta$Site == "Cen16"),]
meta$Site <- as.character(meta$Site)

#read in distance matrix
arsM=read.delim(file = paste(wd, "/data/arsM_rformat_dist_0.03.txt", sep=""))

#add row names back
rownames(arsM)=arsM[,1]

#remove first column
arsM=arsM[,-1]

#make data matrix
arsM=data.matrix(arsM)

#remove first column
arsM=arsM[,-1]


#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#otu table
otu.arsM=otu_table(arsM, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu.arsM, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = FALSE)

#rarefy
rare.arsM=rarefy_even_depth(otu.arsM, sample.size = min(sample_sums(otu.arsM)), rngseed = TRUE)

#check curve
rarecurve(rare.arsM, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.arsM=merge_phyloseq(rare.arsM, metad)

#plot phylo richness
(richness.arsM=plot_richness(phylo.arsM, x="Classification", 
                             color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(richness.arsM, filename = paste(wd, "/figures/arsM.richness.png", sep=""), 
       width = 15)

#calculate evenness
s.arsM=specnumber(rare.arsM)
h.arsM=diversity(rare.arsM, index="shannon")
plieou.arsM=h.arsM/log(s.arsM)

#save evenness number
write.table(plieou.arsM, file = paste(wd, "/output/arsM.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou.arsM=data.frame(plieou.arsM)

#add site column to evenness data
plieou.arsM$Site=rownames(plieou.arsM)

#merge evenness information with fire classification
plieou.arsM=inner_join(plieou.arsM, meta)

#plot evenness by fire classification
(evenness.arsM <- ggplot(plieou.arsM, aes(x = Classification, y = plieou.arsM)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3, width = 0.15) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(evenness.arsM, filename = paste(wd, "/figures/arsM.evenness.png", sep=""), 
       width = 7, height = 5)

#make an output of total gene count per site
arsM.gcounts=rowSums(arsM)

#make object phylo with tree and biom info
tree.arsM <- read.tree(file = paste(wd, "/data/arsM_0.03_tree.nwk", sep=""))
tree.arsM <- phy_tree(tree.arsM)

#merge
phylo.arsM=merge_phyloseq(tree.arsM, rare.arsM, metad)

#plot tree
(arsM.tree.plot <- plot_tree(phylo.arsM, color = "SoilTemperature_to10cm", 
                             size = "abundance",
                             shape = "Classification", label.tips=NULL, 
                             text.size=2, ladderize="left", base.spacing = 0.03, sizebase = 2) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(arsM.tree.plot, filename = paste(wd, "/figures/arsM.tree.png", sep=""), height = 10)

#plot Bray Curtis ordination
ord.arsM <- ordinate(phylo.arsM, method="PCoA", distance="bray")
(bc.ord.arsM=plot_ordination(phylo.arsM, ord.arsM, shape="Classification",
                             title="Bray Curtis") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord.arsM, filename = paste(wd, "/figures/arsM.braycurtis.ord.png", sep=""), 
       width = 6, height = 5)

#plot Sorenson ordination
s.ord.arsM <- ordinate(phylo.arsM, method="PCoA", distance="bray", binary=TRUE)
(sorenson.ord.arsM=plot_ordination(phylo.arsM, s.ord.arsM, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord.arsM, filename = paste(wd, "/figures/arsM.sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#plot unweighted Unifrac ordination
uni.u.ord.arsM <- ordinate(phylo.arsM, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.arsM=plot_ordination(phylo.arsM, uni.u.ord.arsM, shape="Classification", 
                                title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.arsM, filename = paste(wd, "/figures/arsM.u.unifrac.ord.png", 
                                        sep=""), width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.arsM <- ordinate(phylo.arsM, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.arsM=plot_ordination(phylo.arsM, uni.w.ord.arsM, 
                                shape="Classification", title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.arsM, filename = paste(wd, "/figures/arsM.w.unifrac.ord.png", 
                                        sep=""), width = 6, height = 5)


##############################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR arsM#
##############################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names.arsM=list.files(pattern="*arsM_45_taxonabund.txt")
data.arsM <- do.call(rbind, lapply(names.arsM, function(X) {
  data.frame(id = basename(X), read_table(X))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#remove NA rows 
data.arsM <- data.arsM[!is.na(data.arsM$Taxon.Abundance.Fraction.Abundance),]

#split columns 
data.arsM <- data.arsM %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = Taxon.Abundance.Fraction.Abundance, 
           into = c("Taxon", "Abundance", "Fraction.Abundance"), 
           sep = "\t") %>%
  select(-junk) %>%
  group_by(Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.arsM$Fraction.Abundance <- as.numeric(data.arsM$Fraction.Abundance)
data.arsM$Abundance <- as.numeric(data.arsM$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.arsM <- data.arsM %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#change "cen" to "Cen" so it matches outside data
data.arsM$Site <- gsub("cen", "Cen", data.arsM$Site)

#read in arsenic and temperature data
meta <- data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#make column for organism name and join with microbe census data and normalize to it
data.arsM <- data.arsM %>%
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

#get taxonomy for organisms
ncbi <- tax_name(query = data.arsM$organism, get = c("genus", "class", "phylum"), db = "ncbi")

#change query column to "organism"
ncbi$organism <- ncbi$query

#join taxanomic information with abundance information
data.arsM.t <- data.arsM %>%
  left_join(ncbi, by = "organism") %>%
  unique()

#change NA phylum to metagenome
data.arsM.t$phylum[is.na(data.arsM.t$phylum)] = "Metagenome"


#prep colors
n <- 14
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color <- print(sample(col_vector, n))

#order based on temperature
data.arsM.t$Site <- factor(data.arsM.t$Site, 
                           levels = data.arsM.t$Site[order(data.arsM.t$Temp)])
(arsM.abundance.tax.bar.plot <- ggplot(data.arsM.t, aes(x = Site, y = Normalized.Abundance.census, fill = phylum)) +
    geom_bar(stat = "identity") +
    ylab("acr3 Abundance (normalized to genome equivalents)") +
    scale_fill_manual(values = color) +
    theme_classic(base_size = 12))
ggsave(arsM.abundance.tax.bar.plot, filename = paste(wd, "/figures/arsM.phylum.eps", sep=""), width = 10)






#plot data
(arsM.abundance.plot <- ggplot(data.arsM, aes(x = Site, y = Normalized.Abundance.rplB, 
                                              fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow", "green")) +
    ylab("arsM Abundance (normalized to rplB)") +
    theme_classic())

ggsave(arsM.abundance.plot, filename = paste(wd, "/figures/arsM.abundance.png", sep=""))

(arsM.abundance.census.plot <- ggplot(data.arsM, aes(x = Site, 
                                                     y = Normalized.Abundance.census, 
                                                     fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow", "green")) +
    ylab("arsM Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(arsM.abundance.census.plot, filename = paste(wd, "/figures/arsM.abundance.census.png",
                                                    sep=""))


(arsM.abundance.census.taxon.plot <- ggplot(data.arsM, aes(x = organism, 
                                                           y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("arsM abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(arsM.abundance.census.taxon.plot, 
       filename = paste(wd, "/figures/arsM.abundance.census.taxon.png", sep=""), height = 10)

(arsM.abundance.taxon.plot <- ggplot(data.arsM, aes(x = organism, 
                                                    y = Normalized.Abundance.rplB)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("arsM abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(arsM.abundance.taxon.plot, 
       filename = paste(wd, "/figures/arsM.abundance.taxon.png", sep=""), height = 10)

##############################
#arsC_thio DIVERSITY ANALYSIS#
##############################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#remove Cen16 from metadata since we don't have arsC_thio info yet
meta=meta[-which(meta$Site == "Cen01" | meta$Site == "Cen05" | meta$Site == "Cen15" | meta$Site == "Cen16" | meta$Site == "Cen17"),]
meta$Site <- as.character(meta$Site)

#read in distance matrix
arsC_thio=read.delim(file = paste(wd, "/data/arsC_thio_rformat_dist_0.03.txt", sep=""))

#remove Cen15 (too few samples)
arsC_thio <- arsC_thio[-which(arsC_thio$X == "Cen15" | arsC_thio$X == "Cen12"),]

#add row names back
rownames(arsC_thio)=arsC_thio[,1]

#remove first column
arsC_thio=arsC_thio[,-1]

#make data matrix
arsC_thio=data.matrix(arsC_thio)

#remove first column
arsC_thio=arsC_thio[,-1]


#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#otu table
otu.arsC_thio=otu_table(arsC_thio, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu.arsC_thio, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

#rarefy (skip for now)
rare.arsC_thio <- otu.arsC_thio
rare.arsC_thio=rarefy_even_depth(otu.arsC_thio, sample.size = min(sample_sums(otu.arsC_thio)), rngseed = TRUE)

#check curve
rarecurve(rare.arsC_thio, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.arsC_thio=merge_phyloseq(rare.arsC_thio, metad)

#plot phylo richness
(richness.arsC_thio=plot_richness(phylo.arsC_thio, x="Classification", 
                                  color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(richness.arsC_thio, filename = paste(wd, "/figures/arsC_thio.richness.png", sep=""), 
       width = 15)

#calculate evenness
s.arsC_thio=specnumber(rare.arsC_thio)
h.arsC_thio=diversity(rare.arsC_thio, index="shannon")
plieou.arsC_thio=h.arsC_thio/log(s.arsC_thio)

#save evenness number
write.table(plieou.arsC_thio, file = paste(wd, "/output/arsC_thio.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou.arsC_thio=data.frame(plieou.arsC_thio)

#add site column to evenness data
plieou.arsC_thio$Site=rownames(plieou.arsC_thio)

#merge evenness information with fire classification
plieou.arsC_thio=inner_join(plieou.arsC_thio, meta)

#plot evenness by fire classification
(evenness.arsC_thio <- ggplot(plieou.arsC_thio, aes(x = Classification, y = plieou.arsC_thio)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3, width = 0.15) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(evenness.arsC_thio, filename = paste(wd, "/figures/arsC_thio.evenness.png", sep=""), 
       width = 7, height = 5)

#make an output of total gene count per site
arsC_thio.gcounts=rowSums(arsC_thio)

#make object phylo with tree and biom info
tree.arsC_thio <- read.tree(file = paste(wd, "/data/arsC_thio_0.03_tree.nwk", sep=""))
tree.arsC_thio <- phy_tree(tree.arsC_thio)

#merge
phylo.arsC_thio=merge_phyloseq(tree.arsC_thio, rare.arsC_thio, metad)

#plot tree
(arsC_thio.tree.plot <- plot_tree(phylo.arsC_thio, color = "SoilTemperature_to10cm", 
                                  size = "abundance",
                                  shape = "Classification", label.tips="taxa_names", 
                                  text.size=2, ladderize="left", base.spacing = 0.03, sizebase = 2) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(arsC_thio.tree.plot, filename = paste(wd, "/figures/arsC_thio.tree.png", sep=""))

#plot Bray Curtis ordination
ord.arsC_thio <- ordinate(phylo.arsC_thio, method="PCoA", distance="bray")
(bc.ord.arsC_thio=plot_ordination(phylo.arsC_thio, ord.arsC_thio, shape="Classification",
                                  title="Bray Curtis") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord.arsC_thio, filename = paste(wd, "/figures/arsC_thio.braycurtis.ord.png", sep=""), 
       width = 6, height = 5)

#plot Sorenson ordination
s.ord.arsC_thio <- ordinate(phylo.arsC_thio, method="PCoA", distance="bray", binary=TRUE)
(sorenson.ord.arsC_thio=plot_ordination(phylo.arsC_thio, s.ord.arsC_thio, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord.arsC_thio, filename = paste(wd, "/figures/arsC_thio.sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#plot unweighted Unifrac ordination
uni.u.ord.arsC_thio <- ordinate(phylo.arsC_thio, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.arsC_thio=plot_ordination(phylo.arsC_thio, uni.u.ord.arsC_thio, shape="Classification", 
                                     title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.arsC_thio, filename = paste(wd, "/figures/arsC_thio.u.unifrac.ord.png", 
                                             sep=""), width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.arsC_thio <- ordinate(phylo.arsC_thio, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.arsC_thio=plot_ordination(phylo.arsC_thio, uni.w.ord.arsC_thio, 
                                     shape="Classification", title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.arsC_thio, filename = paste(wd, "/figures/arsC_thio.w.unifrac.ord.png", 
                                             sep=""), width = 6, height = 5)


###################################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR arsC_thio#
###################################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names.arsC_thio=list.files(pattern="*arsC_thio_45_taxonabund.txt")
data.arsC_thio <- do.call(rbind, lapply(names.arsC_thio, function(X) {
  data.frame(id = basename(X), read_table(X))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#remove NA rows 
data.arsC_thio <- data.arsC_thio[!is.na(data.arsC_thio$Taxon.Abundance.Fraction.Abundance),]

#split columns 
data.arsC_thio <- data.arsC_thio %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = Taxon.Abundance.Fraction.Abundance, 
           into = c("Taxon", "Abundance", "Fraction.Abundance"), 
           sep = "\t") %>%
  select(-junk) %>%
  group_by(Site)
#for now remove cen12
data.arsC_thio <- data.arsC_thio[-which(data.arsC_thio$Site == "cen12"),]

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.arsC_thio$Fraction.Abundance <- as.numeric(data.arsC_thio$Fraction.Abundance)
data.arsC_thio$Abundance <- as.numeric(data.arsC_thio$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.arsC_thio <- data.arsC_thio %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#change "cen" to "Cen" so it matches outside data
data.arsC_thio$Site <- gsub("cen", "Cen", data.arsC_thio$Site)

#read in arsenic and temperature data
meta <- data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#make column for organism name and join with microbe census data and normalize to it
data.arsC_thio <- data.arsC_thio %>%
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

#plot data
(arsC_thio.abundance.plot <- ggplot(data.arsC_thio, aes(x = Site, y = Normalized.Abundance.rplB, 
                                                        fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow", "green")) +
    ylab("arsC_thio Abundance (normalized to rplB)") +
    theme_classic())

ggsave(arsC_thio.abundance.plot, filename = paste(wd, "/figures/arsC_thio.abundance.png", sep=""))

(arsC_thio.abundance.census.plot <- ggplot(data.arsC_thio, aes(x = Site, 
                                                               y = Normalized.Abundance.census, 
                                                               fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow")) +
    ylab("arsC_thio Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(arsC_thio.abundance.census.plot, filename = paste(wd, "/figures/arsC_thio.abundance.census.png",
                                                         sep=""))


(arsC_thio.abundance.census.taxon.plot <- ggplot(data.arsC_thio, aes(x = organism, 
                                                                     y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("arsC_thio abundance (normalized to genome equivalents)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(arsC_thio.abundance.census.taxon.plot, 
       filename = paste(wd, "/figures/arsC_thio.abundance.census.taxon.png", sep=""), height = 10)

(arsC_thio.abundance.taxon.plot <- ggplot(data.arsC_thio, aes(x = organism, 
                                                              y = Normalized.Abundance.rplB)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("arsC_thio abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(arsC_thio.abundance.taxon.plot, 
       filename = paste(wd, "/figures/arsC_thio.abundance.taxon.png", sep=""), height = 10)

#summarise arsC_thio data
data.arsC_thio.sum <- data.arsC_thio %>%
  group_by(Classification, Site) %>%
  summarise(Abundance.census.total = sum(Normalized.Abundance.census), 
            Abundance.rplB.total = sum(Normalized.Abundance.rplB))

(arsC_thio.abundance.census.box <- ggplot(data.arsC_thio.sum, 
                                          aes(x = Classification, 
                                              y = Abundance.census.total)) +
  geom_boxplot() +
    geom_jitter())

##############################
#arsC_glut DIVERSITY ANALYSIS#
##############################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#remove Cen16 from metadata since we don't have arsC_glut info yet
meta=meta[which(meta$Site == "Cen07" | meta$Site == "Cen10" | meta$Site == "Cen01" | meta$Site == "Cen03" | meta$Site == "Cen12" | meta$Site == "Cen04" | meta$Site == "Cen06"),]
meta$Site <- as.character(meta$Site)

#read in distance matrix
arsC_glut=read.delim(file = paste(wd, "/data/arsC_glut_rformat_dist_0.03.txt", sep=""))

#add row names back
rownames(arsC_glut)=arsC_glut[,1]

#remove first column
arsC_glut=arsC_glut[,-1]

#make data matrix
arsC_glut=data.matrix(arsC_glut)

#remove first column
arsC_glut=arsC_glut[,-1]


#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#otu table
otu.arsC_glut=otu_table(arsC_glut, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu.arsC_glut, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = FALSE)

#rarefy
rare.arsC_glut=rarefy_even_depth(otu.arsC_glut, sample.size = min(sample_sums(otu.arsC_glut)), rngseed = TRUE)

#check curve
rarecurve(rare.arsC_glut, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.arsC_glut=merge_phyloseq(rare.arsC_glut, metad)

#plot phylo richness
(richness.arsC_glut=plot_richness(phylo.arsC_glut, x="Classification", 
                                  color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(richness.arsC_glut, filename = paste(wd, "/figures/arsC_glut.richness.png", sep=""), 
       width = 15)

#calculate evenness
s.arsC_glut=specnumber(rare.arsC_glut)
h.arsC_glut=diversity(rare.arsC_glut, index="shannon")
plieou.arsC_glut=h.arsC_glut/log(s.arsC_glut)

#save evenness number
write.table(plieou.arsC_glut, file = paste(wd, "/output/arsC_glut.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou.arsC_glut=data.frame(plieou.arsC_glut)

#add site column to evenness data
plieou.arsC_glut$Site=rownames(plieou.arsC_glut)

#merge evenness information with fire classification
plieou.arsC_glut=inner_join(plieou.arsC_glut, meta)

#plot evenness by fire classification
(evenness.arsC_glut <- ggplot(plieou.arsC_glut, aes(x = Classification, y = plieou.arsC_glut)) +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3, width = 0.15) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(evenness.arsC_glut, filename = paste(wd, "/figures/arsC_glut.evenness.png", sep=""), 
       width = 7, height = 5)

#make an output of total gene count per site
arsC_glut.gcounts=rowSums(arsC_glut)

#make object phylo with tree and biom info
tree.arsC_glut <- read.tree(file = paste(wd, "/data/arsC_glut_0.03_tree.nwk", sep=""))
tree.arsC_glut <- phy_tree(tree.arsC_glut)

#merge
phylo.arsC_glut=merge_phyloseq(tree.arsC_glut, rare.arsC_glut, metad)

#plot tree
(arsC_glut.tree.plot <- plot_tree(phylo.arsC_glut, color = "SoilTemperature_to10cm", 
                                  size = "abundance",
                                  shape = "Classification", label.tips="taxa_names", 
                                  text.size=2, ladderize="left", base.spacing = 0.03, sizebase = 2) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(arsC_glut.tree.plot, filename = paste(wd, "/figures/arsC_glut.tree.png", sep=""))

#plot Bray Curtis ordination
ord.arsC_glut <- ordinate(phylo.arsC_glut, method="PCoA", distance="bray")
(bc.ord.arsC_glut=plot_ordination(phylo.arsC_glut, ord.arsC_glut, shape="Classification",
                                  title="Bray Curtis") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord.arsC_glut, filename = paste(wd, "/figures/arsC_glut.braycurtis.ord.png", sep=""), 
       width = 6, height = 5)

#plot Sorenson ordination
s.ord.arsC_glut <- ordinate(phylo.arsC_glut, method="PCoA", distance="bray", binary=TRUE)
(sorenson.ord.arsC_glut=plot_ordination(phylo.arsC_glut, s.ord.arsC_glut, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord.arsC_glut, filename = paste(wd, "/figures/arsC_glut.sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#plot unweighted Unifrac ordination
uni.u.ord.arsC_glut <- ordinate(phylo.arsC_glut, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.arsC_glut=plot_ordination(phylo.arsC_glut, uni.u.ord.arsC_glut, shape="Classification", 
                                     title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.arsC_glut, filename = paste(wd, "/figures/arsC_glut.u.unifrac.ord.png", 
                                             sep=""), width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.arsC_glut <- ordinate(phylo.arsC_glut, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.arsC_glut=plot_ordination(phylo.arsC_glut, uni.w.ord.arsC_glut, 
                                     shape="Classification", title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.arsC_glut, filename = paste(wd, "/figures/arsC_glut.w.unifrac.ord.png", 
                                             sep=""), width = 6, height = 5)


###################################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR arsC_glut#
###################################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names.arsC_glut=list.files(pattern="*arsC_glut_45_taxonabund.txt")
data.arsC_glut <- do.call(rbind, lapply(names.arsC_glut, function(X) {
  data.frame(id = basename(X), read_table(X))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#remove NA rows 
data.arsC_glut <- data.arsC_glut[!is.na(data.arsC_glut$Taxon.Abundance.Fraction.Abundance),]

#split columns 
data.arsC_glut <- data.arsC_glut %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = Taxon.Abundance.Fraction.Abundance, 
           into = c("Taxon", "Abundance", "Fraction.Abundance"), 
           sep = "\t") %>%
  select(-junk) %>%
  group_by(Site)
#for now remove cen12
data.arsC_glut <- data.arsC_glut[-which(data.arsC_glut$Site == "cen12"),]

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.arsC_glut$Fraction.Abundance <- as.numeric(data.arsC_glut$Fraction.Abundance)
data.arsC_glut$Abundance <- as.numeric(data.arsC_glut$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.arsC_glut <- data.arsC_glut %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#change "cen" to "Cen" so it matches outside data
data.arsC_glut$Site <- gsub("cen", "Cen", data.arsC_glut$Site)

#read in arsenic and temperature data
meta <- data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#make column for organism name and join with microbe census data and normalize to it
data.arsC_glut <- data.arsC_glut %>%
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
data.arsC_glut <- data.arsC_glut[-which(data.arsC_glut$Site == "Cen01"),]

#plot data
(arsC_glut.abundance.plot <- ggplot(data.arsC_glut, aes(x = Site, y = Normalized.Abundance.rplB, 
                                                        fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow")) +
    ylab("arsC_glut Abundance (normalized to rplB)") +
    theme_classic())

ggsave(arsC_glut.abundance.plot, filename = paste(wd, "/figures/arsC_glut.abundance.png", sep=""))

(arsC_glut.abundance.census.plot <- ggplot(data.arsC_glut, aes(x = Site, 
                                                               y = Normalized.Abundance.census, 
                                                               fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "yellow")) +
    ylab("arsC_glut Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(arsC_glut.abundance.census.plot, filename = paste(wd, "/figures/arsC_glut.abundance.census.png",
                                                         sep=""))


(arsC_glut.abundance.census.taxon.plot <- ggplot(data.arsC_glut, aes(x = organism, 
                                                                     y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("arsC_glut abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(arsC_glut.abundance.census.taxon.plot, 
       filename = paste(wd, "/figures/arsC_glut.abundance.census.taxon.png", sep=""), height = 10)

(arsC_glut.abundance.taxon.plot <- ggplot(data.arsC_glut, aes(x = organism, 
                                                              y = Normalized.Abundance.rplB)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("arsC_glut abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(arsC_glut.abundance.taxon.plot, 
       filename = paste(wd, "/figures/arsC_glut.abundance.taxon.png", sep=""), height = 10)


