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
#ClassC DIVERSITY ANALYSIS#
#########################

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#read in distance matrix
ClassC=read.delim(file = paste(wd, "/data/ClassC_rformat_dist_0.03.txt", sep=""))

#add row names back
rownames(ClassC)=ClassC[,1]

#remove first column
ClassC=ClassC[,-1]

#make data matrix
ClassC=data.matrix(ClassC)

#remove first column
ClassC=ClassC[,-1]

#call metadata sample data
metad=meta[-1]
rownames(metad)=meta$Site
metad=sample_data(metad)

#otu table
otu.ClassC=otu_table(ClassC, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu.ClassC, step=5, col = rarecol, label = TRUE)

#rarefy
rare.ClassC=rarefy_even_depth(otu.ClassC, sample.size = min(sample_sums(otu.ClassC)), rngseed = TRUE)

#check curve
rarecurve(rare.ClassC, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.ClassC=merge_phyloseq(rare.ClassC, metad)

#plot phylo richness
(ClassC.richness=plot_richness(phylo.ClassC, x="Classification", color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(ClassC.richness, filename = paste(wd, "/figures/ClassC.richness.png", sep=""), 
       width = 15)

#calculate evenness
s.ClassC=specnumber(rare.ClassC)
h.ClassC=diversity(rare.ClassC, index="shannon")
plieou.ClassC=h.ClassC/log(s.ClassC)

#save evenness number
write.table(plieou.ClassC, file = paste(wd, "/output/ClassC.evenness.txt", sep=""))

#make plieou a dataframe for plotting
plieou.ClassC=data.frame(plieou.ClassC)

#add site column to evenness data
plieou.ClassC$Site=rownames(plieou.ClassC)

#merge evenness information with fire classification
plieou.ClassC=inner_join(plieou.ClassC, meta)

#plot evenness by fire classification
(ClassC.evenness <- ggplot(plieou.ClassC, aes(x = Classification, y = plieou.ClassC)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(ClassC.evenness, filename = paste(wd, "/figures/ClassC.evenness.png", sep=""), 
       width = 7, height = 5)

#make an output of total gene count per site
ClassC.gcounts=rowSums(ClassC)

#make object phylo with tree and biom info
tree.ClassC <- read.tree(file = paste(wd, "/data/ClassC_0.03_tree.nwk", sep=""))
tree.ClassC <- phy_tree(tree.ClassC)

#merge
phylo.ClassC=merge_phyloseq(tree.ClassC, rare.ClassC, metad)


#plot Bray Curtis ordination
ord.ClassC <- ordinate(phylo.ClassC, method="PCoA", distance="bray")
(bc.ord.ClassC=plot_ordination(phylo.ClassC, ord.ClassC, shape="Classification", title="Bray Curtis") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord.ClassC, filename = paste(wd, "/figures/ClassC.braycurtis.ord.png", sep=""), 
       width = 6, height = 5)

#plot Sorenson ordination
s.ord.ClassC <- ordinate(phylo.ClassC, method="PCoA", distance="bray", binary=TRUE)
(sorenson.ord.ClassC=plot_ordination(phylo.ClassC, s.ord.ClassC, shape="Classification", title="Sorenson") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save sorenson ordination
ggsave(sorenson.ord.ClassC, filename = paste(wd, "/figures/ClassC.sorenson.ord.png", sep=""), 
       width = 6, height = 5)

#plot unweighted Unifrac ordination
uni.u.ord.ClassC <- ordinate(phylo.ClassC, method="PCoA", distance="unifrac", weighted = FALSE)
(uni.u.ord.ClassC=plot_ordination(phylo.ClassC, uni.u.ord.ClassC, shape="Classification", 
                                title="Unweighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save unweighted Unifrac ordination
ggsave(uni.u.ord.ClassC, filename = paste(wd, "/figures/ClassC.u.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot weighted Unifrac ordination
uni.w.ord.ClassC <- ordinate(phylo.ClassC, method="PCoA", distance="unifrac", weighted = TRUE)
(uni.w.ord.ClassC=plot_ordination(phylo.ClassC, uni.w.ord.ClassC, shape="Classification", 
                                title="Weighted Unifrac") +
    geom_point(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save weighted Unifrac ordination
ggsave(uni.w.ord.ClassC, filename = paste(wd, "/figures/ClassC.w.unifrac.ord.png", sep=""), 
       width = 6, height = 5)

#plot tree
(ClassC.tree.plot <- plot_tree(phylo.ClassC, color = "SoilTemperature_to10cm", 
                             size = "abundance",
                             shape = "Classification", label.tips=NULL, 
                             text.size=2, ladderize="left", base.spacing = 0.04) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(ClassC.tree.plot, filename = paste(wd, "/figures/ClassC.tree.png", sep=""), 
       height = 10, width = 10)

##############################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR ClassC#
##############################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names.ClassC=list.files(pattern="*ClassC_45_taxonabund.txt")
data.ClassC <- do.call(rbind, lapply(names.ClassC, function(X) {
  data.frame(id = basename(X), read_table(X))}))

#move back up a directory to proceed with analysis
setwd("../")
wd <- print(getwd())

#remove NA rows 
data.ClassC <- data.ClassC[!is.na(data.ClassC$Taxon.Abundance.Fraction.Abundance),]

#split columns 
data.ClassC <- data.ClassC %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = Taxon.Abundance.Fraction.Abundance, 
           into = c("Taxon", "Abundance", "Fraction.Abundance"), 
           sep = "\t") %>%
  select(-junk) %>%
  group_by(Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.ClassC$Fraction.Abundance <- as.numeric(data.ClassC$Fraction.Abundance)
data.ClassC$Abundance <- as.numeric(data.ClassC$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.ClassC <- data.ClassC %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#change "cen" to "Cen" so it matches outside data
data.ClassC$Site <- gsub("cen", "Cen", data.ClassC$Site)

#read in arsenic and temperature data
meta <- data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#make column for organism name and join with microbe census data and normalize to it
data.ClassC <- data.ClassC %>%
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
ClassC.ncbi <- tax_name(query = data.ClassC$organism, get = c("genus", "class", "phylum"), db = "ncbi")

#change query column to "organism"
ClassC.ncbi$organism <- ClassC.ncbi$query

#join taxanomic information with abundance information
data.ClassC.t <- data.ClassC %>%
  left_join(ClassC.ncbi, by = "organism") %>%
  unique()

#change NA phylum to metagenome
data.ClassC.t$phylum[is.na(data.ClassC.t$phylum)] = "Metagenome"

#order based on temperature
data.ClassC.t$Site <- factor(data.ClassC.t$Site, 
                           levels = data.ClassC.t$Site[order(data.ClassC.t$Temp)])

#prep colors
color.ClassC <- c("#FF7F00", "#7570B3", "#CAB2D6", "#FBB4AE", "#F0027F", "#BEBADA", "#E78AC3", "#A6D854", "#B3B3B3", "#386CB0", "#BC80BD", "#FFFFCC", "#BF5B17", "#984EA3", "#CCCCCC", "#FFFF99", "#B15928", "#F781BF", "#FDC086", "#A6CEE3", "#FDB462", "#FED9A6", "#E6AB02", "#E31A1C", "#B2DF8A", "#377EB8", "#FCCDE5", "#80B1D3", "#FFD92F", "#33A02C", "#66C2A5", "#666666")


(ClassC.abundance.tax.bar.plot <- ggplot(data.ClassC.t, aes(x = Site, y = Normalized.Abundance.census, fill = phylum)) +
    geom_bar(stat = "identity") +
    ylab("ClassC Abundance (normalized to genome equivalents)") +
    scale_fill_manual(values = color.ClassC) +
    theme_classic(base_size = 12))
ggsave(ClassC.abundance.tax.bar.plot, filename = paste(wd, "/figures/ClassC.census.phylum.eps", sep=""), width = 10)

#order based on temperature
data.ClassC$Site <- factor(data.ClassC$Site, 
                         levels = data.ClassC$Site[order(data.ClassC$Temp)])

#plot data
(ClassC.abundance.plot <- ggplot(data.ClassC.t, aes(x = Site, y = Normalized.Abundance.rplB, 
                                              fill = phylum)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color.ClassC) +
    ylab("ClassC Abundance (normalized to rplB)") +
    theme_classic())

ggsave(ClassC.abundance.plot, filename = paste(wd, "/figures/ClassC.rplB.phylum.png", sep=""))

(ClassC.abundance.census.plot <- ggplot(data.ClassC, aes(x = Site, 
                                                     y = Normalized.Abundance.census, 
                                                     fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("firebrick2", "yellow1", "green3")) +
    ylab("ClassC Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(ClassC.abundance.census.plot, 
       filename = paste(wd, "/figures/ClassC.abundance.census.png", sep=""))


(ClassC.abundance.census.taxon.plot <- ggplot(data.ClassC, aes(x = organism, 
                                                           y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("ClassC abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(ClassC.abundance.census.taxon.plot, 
       filename = paste(wd, "/figures/ClassC.abundance.census.taxon.png", sep=""), height = 10)

(ClassC.abundance.taxon.plot <- ggplot(data.ClassC, aes(x = organism, 
                                                    y = Normalized.Abundance.rplB)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("ClassC abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))
ggsave(ClassC.abundance.taxon.plot, 
       filename = paste(wd, "/figures/ClassC.abundance.taxon.png", sep=""), height = 10)

#summarise ClassC information
data.ClassC.sum <- data.ClassC %>%
  group_by(Classification, Site, Temp) %>%
  summarise(Abundance.census.tot = sum(Normalized.Abundance.census), 
            Abundance.rplB.tot = sum(Normalized.Abundance.rplB))

(ClassC.abundance.census.box <- ggplot(data.ClassC.sum, aes(x = Classification, 
                                                        y = Abundance.census.tot)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, aes(color = Temp)) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
  ylab("ClassC abundance (normalized to genome equivalents)"))

ggsave(ClassC.abundance.census.box, 
       filename = paste(wd, "/figures/ClassC.abundance.census.box.png", sep=""))

(ClassC.abundance.rplB.box <- ggplot(data.ClassC.sum, aes(x = Classification, 
                                                        y = Abundance.rplB.tot)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, aes(color = Temp)) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("ClassC abundance (normalized to genome equivalents)"))

ggsave(ClassC.abundance.rplB.box, 
       filename = paste(wd, "/figures/ClassC.abundance.rplB.box.png", sep=""))

