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
#vanH DIVERSITY ANALYSIS#
#########################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#read in distance matrix
vanH=read.delim(file = paste(wd, "/data/vanH_rformat_dist_0.03.txt", sep=""))

#add row names back
rownames(vanH)=vanH[,1]

#remove first column
vanH=vanH[,-1]

#make data matrix
vanH=data.matrix(vanH)

#remove first column
vanH=vanH[,-1]

#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#otu table
otu.vanH=otu_table(vanH, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu.vanH, step=5, col = rarecol, label = TRUE)

#rarefy
rare.vanH=rarefy_even_depth(otu.vanH, sample.size = min(sample_sums(otu.vanH)), rngseed = TRUE)

#check curve
rarecurve(rare.vanH, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.vanH=merge_phyloseq(rare.vanH, metad)

#plot phylo richness
(vanH.richness=plot_richness(phylo.vanH, x="Classification", color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(vanH.richness, filename = paste(wd, "/figures/vanH.richness.png", sep=""), 
       width = 15)

#calculate evenness
s.vanH=specnumber(rare.vanH)
h.vanH=diversity(rare.vanH, index="shannon")
plieou.vanH=h.vanH/log(s.vanH)

#save evenness number
write.table(plieou.vanH, file = paste(wd, "/output/vanH.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou.vanH=data.frame(plieou.vanH)

#add site column to evenness data
plieou.vanH$Site=rownames(plieou.vanH)

#merge evenness information with fire classification
plieou.vanH=inner_join(plieou.vanH, meta)

#plot evenness by fire classification
(vanH.evenness <- ggplot(plieou.vanH, aes(x = Classification, y = plieou.vanH)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(vanH.evenness, filename = paste(wd, "/figures/vanH.evenness.png", sep=""), 
       width = 7, height = 5)

#make an output of total gene count per site
vanH.gcounts=rowSums(vanH)

#make object phylo with tree and biom info
tree.vanH <- read.tree(file = paste(wd, "/data/vanH_0.03_tree.nwk", sep=""))
tree.vanH <- phy_tree(tree.vanH)

#merge
phylo.vanH=merge_phyloseq(tree.vanH, rare.vanH, metad)


#plot Bray Curtis ordination
ord.vanH <- ordinate(phylo.vanH, method="PCoA", distance="bray")
(bc.ord.vanH=plot_ordination(phylo.vanH, ord.vanH, shape="Classification", title="Bray Curtis") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord.vanH, filename = paste(wd, "/figures/vanH.braycurtis.ord.png", sep=""), 
       width = 6, height = 5)

#plot Sorenson ordination
s.ord.vanH <- ordinate(phylo.vanH, method="PCoA", distance="bray", binary=TRUE)
(sorenson.ord.vanH=plot_ordination(phylo.vanH, s.ord.vanH, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord.vanH, filename = paste(wd, "/figures/vanH.sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#plot unweighted Unifrac ordination
uni.u.ord.vanH <- ordinate(phylo.vanH, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.vanH=plot_ordination(phylo.vanH, uni.u.ord.vanH, shape="Classification", 
                                title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.vanH, filename = paste(wd, "/figures/vanH.u.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.vanH <- ordinate(phylo.vanH, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.vanH=plot_ordination(phylo.vanH, uni.w.ord.vanH, shape="Classification", 
                                title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.vanH, filename = paste(wd, "/figures/vanH.w.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot tree
(vanH.tree.plot <- plot_tree(phylo.vanH, color = "SoilTemperature_to10cm", 
                             size = "abundance",
                             shape = "Classification", label.tips=NULL, 
                             text.size=2, ladderize="left", base.spacing = 0.04) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(vanH.tree.plot, filename = paste(wd, "/figures/vanH.tree.png", sep=""), 
       height = 10, width = 10)

##############################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR vanH#
##############################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names.vanH=list.files(pattern="*vanH_45_taxonabund.txt")
data.vanH <- do.call(rbind, lapply(names.vanH, function(X) {
  data.frame(id = basename(X), read_table(X))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#remove NA rows 
data.vanH <- data.vanH[!is.na(data.vanH$Taxon.Abundance.Fraction.Abundance),]

#split columns 
data.vanH <- data.vanH %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = Taxon.Abundance.Fraction.Abundance, 
           into = c("Taxon", "Abundance", "Fraction.Abundance"), 
           sep = "\t") %>%
  select(-junk) %>%
  group_by(Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.vanH$Fraction.Abundance <- as.numeric(data.vanH$Fraction.Abundance)
data.vanH$Abundance <- as.numeric(data.vanH$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.vanH <- data.vanH %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#change "cen" to "Cen" so it matches outside data
data.vanH$Site <- gsub("cen", "Cen", data.vanH$Site)

#read in arsenic and temperature data
meta <- data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#make column for organism name and join with microbe census data and normalize to it
data.vanH <- data.vanH %>%
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
vanH.ncbi <- tax_name(query = data.vanH$organism, get = c("genus", "class", "phylum"), db = "ncbi")

#change query column to "organism"
vanH.ncbi$organism <- vanH.ncbi$query

#join taxanomic information with abundance information
data.vanH.t <- data.vanH %>%
  left_join(vanH.ncbi, by = "organism") %>%
  unique()

#change NA phylum to metagenome
data.vanH.t$phylum[is.na(data.vanH.t$phylum)] = "Metagenome"

#order based on temperature
data.vanH.t$Site <- factor(data.vanH.t$Site, 
                           levels = data.vanH.t$Site[order(data.vanH.t$Temp)])

#prep colors
color.vanH <- c("#FF7F00", "#7570B3", "#CAB2D6", "#FBB4AE", "#F0027F", "#BEBADA", "#E78AC3", "#A6D854", "#B3B3B3", "#386CB0", "#BC80BD", "#FFFFCC", "#BF5B17", "#984EA3", "#CCCCCC", "#FFFF99", "#B15928", "#F781BF", "#FDC086", "#A6CEE3", "#FDB462", "#FED9A6", "#E6AB02", "#E31A1C", "#B2DF8A", "#377EB8", "#FCCDE5", "#80B1D3", "#FFD92F", "#33A02C", "#66C2A5", "#666666")


(vanH.abundance.tax.bar.plot <- ggplot(data.vanH.t, aes(x = Site, y = Normalized.Abundance.census, fill = phylum)) +
    geom_bar(stat = "identity") +
    ylab("vanH Abundance (normalized to genome equivalents)") +
    scale_fill_manual(values = color.vanH) +
    theme_classic(base_size = 12))
ggsave(vanH.abundance.tax.bar.plot, filename = paste(wd, "/figures/vanH.census.phylum.eps", sep=""), width = 10)

#order based on temperature
data.vanH$Site <- factor(data.vanH$Site, 
                         levels = data.vanH$Site[order(data.vanH$Temp)])

#plot data
(vanH.abundance.plot <- ggplot(data.vanH.t, aes(x = Site, y = Normalized.Abundance.rplB, 
                                              fill = phylum)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color.vanH) +
    ylab("vanH Abundance (normalized to rplB)") +
    theme_classic())

ggsave(vanH.abundance.plot, filename = paste(wd, "/figures/vanH.rplB.phylum.png", sep=""))

(vanH.abundance.census.plot <- ggplot(data.vanH, aes(x = Site, 
                                                     y = Normalized.Abundance.census, 
                                                     fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("firebrick2", "yellow1", "green3")) +
    ylab("vanH Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(vanH.abundance.census.plot, 
       filename = paste(wd, "/figures/vanH.abundance.census.png", sep=""))


(vanH.abundance.census.taxon.plot <- ggplot(data.vanH, aes(x = organism, 
                                                           y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("vanH abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(vanH.abundance.census.taxon.plot, 
       filename = paste(wd, "/figures/vanH.abundance.census.taxon.png", sep=""), height = 10)

(vanH.abundance.taxon.plot <- ggplot(data.vanH, aes(x = organism, 
                                                    y = Normalized.Abundance.rplB)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("vanH abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(vanH.abundance.taxon.plot, 
       filename = paste(wd, "/figures/vanH.abundance.taxon.png", sep=""), height = 10)

#summarise vanH information
data.vanH.sum <- data.vanH %>%
  group_by(Classification, Site, Temp) %>%
  summarise(Abundance.census.tot = sum(Normalized.Abundance.census), 
            Abundance.rplB.tot = sum(Normalized.Abundance.rplB))

(vanH.abundance.census.box <- ggplot(data.vanH.sum, aes(x = Classification, 
                                                        y = Abundance.census.tot)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, aes(color = Temp)) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
  ylab("vanH abundance (normalized to genome equivalents)"))

ggsave(vanH.abundance.census.box, 
       filename = paste(wd, "/figures/vanH.abundance.census.box.png", sep=""))

(vanH.abundance.census.box <- ggplot(data.vanH.sum, aes(x = Classification, 
                                                        y = Abundance.rplB.tot)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, aes(color = Temp)) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("vanH abundance (normalized to genome equivalents)"))

ggsave(vanH.abundance.census.box, 
       filename = paste(wd, "/figures/vanH.abundance.rplB.box.png", sep=""))
