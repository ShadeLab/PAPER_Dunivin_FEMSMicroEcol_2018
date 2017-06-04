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

#make colores for rarefaction curves (n=12)
rarecol <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink", "green", "red", "brown", "grey", "purple")

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
rarecurve(otu, step=5, col = rarecol, label = FALSE)

#rarefy
rare=rarefy_even_depth(otu, sample.size = min(sample_sums(otu)), rngseed = TRUE)

#check curve
rarecurve(rare, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo=merge_phyloseq(rare, metad)

#plot phylo richness
(richness=plot_richness(phylo, x="Classification", color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

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
#intI DIVERSITY ANALYSIS#
#########################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#read in distance matrix
intI=read.delim(file = paste(wd, "/data/intI_rformat_dist_0.03.txt", sep=""))

#add row names back
rownames(intI)=intI[,1]

#remove first column
intI=intI[,-1]

#make data matrix
intI=data.matrix(intI)

#remove first column
intI=intI[,-1]

#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#otu table
otu.intI=otu_table(intI, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu.intI, step=5, col = rarecol, label = TRUE)

#rarefy
rare.intI=rarefy_even_depth(otu.intI, sample.size = min(sample_sums(otu.intI)), rngseed = TRUE)

#check curve
rarecurve(rare.intI, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.intI=merge_phyloseq(rare.intI, metad)

#plot phylo richness
(intI.richness=plot_richness(phylo.intI, x="Classification", color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(intI.richness, filename = paste(wd, "/figures/intI.richness.png", sep=""), 
       width = 15)

#calculate evenness
s.intI=specnumber(rare.intI)
h.intI=diversity(rare.intI, index="shannon")
plieou.intI=h.intI/log(s.intI)

#save evenness number
write.table(plieou.intI, file = paste(wd, "/output/intI.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou.intI=data.frame(plieou.intI)

#add site column to evenness data
plieou.intI$Site=rownames(plieou.intI)

#merge evenness information with fire classification
plieou.intI=inner_join(plieou.intI, meta)

#plot evenness by fire classification
(intI.evenness <- ggplot(plieou.intI, aes(x = Classification, y = plieou.intI)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(intI.evenness, filename = paste(wd, "/figures/intI.evenness.png", sep=""), 
       width = 7, height = 5)

#make an output of total gene count per site
intI.gcounts=rowSums(intI)

#make object phylo with tree and biom info
tree.intI <- read.tree(file = paste(wd, "/data/intI_0.03_tree.nwk", sep=""))
tree.intI <- phy_tree(tree.intI)

#merge
phylo.intI=merge_phyloseq(tree.intI, rare.intI, metad)


#plot Bray Curtis ordination
ord.intI <- ordinate(phylo.intI, method="PCoA", distance="bray")
(bc.ord.intI=plot_ordination(phylo.intI, ord.intI, shape="Classification", title="Bray Curtis") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord.intI, filename = paste(wd, "/figures/intI.braycurtis.ord.png", sep=""), 
       width = 6, height = 5)

#plot Sorenson ordination
s.ord.intI <- ordinate(phylo.intI, method="PCoA", distance="bray", binary=TRUE)
(sorenson.ord.intI=plot_ordination(phylo.intI, s.ord.intI, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord.intI, filename = paste(wd, "/figures/intI.sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#plot unweighted Unifrac ordination
uni.u.ord.intI <- ordinate(phylo.intI, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.intI=plot_ordination(phylo.intI, uni.u.ord.intI, shape="Classification", 
                                title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.intI, filename = paste(wd, "/figures/intI.u.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.intI <- ordinate(phylo.intI, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.intI=plot_ordination(phylo.intI, uni.w.ord.intI, shape="Classification", 
                                title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.intI, filename = paste(wd, "/figures/intI.w.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot tree
(intI.tree.plot <- plot_tree(phylo.intI, color = "SoilTemperature_to10cm", 
                             size = "abundance",
                             shape = "Classification", label.tips=NULL, 
                             text.size=2, ladderize="left", base.spacing = 0.04) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(intI.tree.plot, filename = paste(wd, "/figures/intI.tree.png", sep=""), 
       height = 10, width = 10)

##############################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR intI#
##############################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names.intI=list.files(pattern="*intI_45_taxonabund.txt")
data.intI <- do.call(rbind, lapply(names.intI, function(X) {
  data.frame(id = basename(X), read_table(X))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#remove NA rows 
data.intI <- data.intI[!is.na(data.intI$Taxon.Abundance.Fraction.Abundance),]

#split columns 
data.intI <- data.intI %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = Taxon.Abundance.Fraction.Abundance, 
           into = c("Taxon", "Abundance", "Fraction.Abundance"), 
           sep = "\t") %>%
  select(-junk) %>%
  group_by(Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.intI$Fraction.Abundance <- as.numeric(data.intI$Fraction.Abundance)
data.intI$Abundance <- as.numeric(data.intI$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.intI <- data.intI %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#change "cen" to "Cen" so it matches outside data
data.intI$Site <- gsub("cen", "Cen", data.intI$Site)

#read in arsenic and temperature data
meta <- data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#make column for organism name and join with microbe census data and normalize to it
data.intI <- data.intI %>%
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
intI.ncbi <- tax_name(query = data.intI$organism, get = c("genus", "class", "phylum"), db = "ncbi")

#change query column to "organism"
intI.ncbi$organism <- intI.ncbi$query

#join taxanomic information with abundance information
data.intI.t <- data.intI %>%
  left_join(intI.ncbi, by = "organism") %>%
  unique()

#change NA phylum to metagenome
data.intI.t$phylum[is.na(data.intI.t$phylum)] = "Metagenome"

#order based on temperature
data.intI.t$Site <- factor(data.intI.t$Site, 
                           levels = data.intI.t$Site[order(data.intI.t$Temp)])

#prep colors
color.intI <- c("#FF7F00", "#7570B3", "#CAB2D6", "#FBB4AE", "#F0027F", "#BEBADA", "#E78AC3", "#A6D854", "#B3B3B3", "#386CB0", "#BC80BD", "#FFFFCC", "#BF5B17", "#984EA3", "#CCCCCC", "#FFFF99", "#B15928", "#F781BF", "#FDC086", "#A6CEE3", "#FDB462", "#FED9A6", "#E6AB02", "#E31A1C", "#B2DF8A", "#377EB8", "#FCCDE5", "#80B1D3", "#FFD92F", "#33A02C", "#66C2A5", "#666666")


(intI.abundance.tax.bar.plot <- ggplot(data.intI.t, aes(x = Site, y = Normalized.Abundance.census, fill = phylum)) +
    geom_bar(stat = "identity") +
    ylab("intI Abundance (normalized to genome equivalents)") +
    scale_fill_manual(values = color.intI) +
    theme_classic(base_size = 12))
ggsave(intI.abundance.tax.bar.plot, filename = paste(wd, "/figures/intI.census.phylum.eps", sep=""), width = 10)

#order based on temperature
data.intI$Site <- factor(data.intI$Site, 
                         levels = data.intI$Site[order(data.intI$Temp)])

#plot data
(intI.abundance.plot <- ggplot(data.intI.t, aes(x = Site, y = Normalized.Abundance.rplB, 
                                              fill = phylum)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color.intI) +
    ylab("intI Abundance (normalized to rplB)") +
    theme_classic())

ggsave(intI.abundance.plot, filename = paste(wd, "/figures/intI.rplB.phylum.png", sep=""))

(intI.abundance.census.plot <- ggplot(data.intI, aes(x = Site, 
                                                     y = Normalized.Abundance.census, 
                                                     fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("firebrick2", "yellow1", "green3")) +
    ylab("intI Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(intI.abundance.census.plot, 
       filename = paste(wd, "/figures/intI.abundance.census.png", sep=""))


(intI.abundance.census.taxon.plot <- ggplot(data.intI, aes(x = organism, 
                                                           y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("intI abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(intI.abundance.census.taxon.plot, 
       filename = paste(wd, "/figures/intI.abundance.census.taxon.png", sep=""), height = 10)

(intI.abundance.taxon.plot <- ggplot(data.intI, aes(x = organism, 
                                                    y = Normalized.Abundance.rplB)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("intI abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(intI.abundance.taxon.plot, 
       filename = paste(wd, "/figures/intI.abundance.taxon.png", sep=""), height = 10)

#summarise intI information
data.intI.sum <- data.intI %>%
  group_by(Classification, Site, Temp) %>%
  summarise(Abundance.census.tot = sum(Normalized.Abundance.census), 
            Abundance.rplB.tot = sum(Normalized.Abundance.rplB))

(intI.abundance.census.box <- ggplot(data.intI.sum, aes(x = Classification, 
                                                        y = Abundance.census.tot)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, aes(color = Temp)) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
  ylab("intI abundance (normalized to genome equivalents)"))

ggsave(intI.abundance.census.box, 
       filename = paste(wd, "/figures/intI.abundance.census.box.png", sep=""))

(intI.abundance.rplB.box <- ggplot(data.intI.sum, aes(x = Classification, 
                                                        y = Abundance.rplB.tot)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, aes(color = Temp)) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("intI abundance (normalized to genome equivalents)"))

ggsave(intI.abundance.rplB.box, 
       filename = paste(wd, "/figures/intI.abundance.rplB.box.png", sep=""))
