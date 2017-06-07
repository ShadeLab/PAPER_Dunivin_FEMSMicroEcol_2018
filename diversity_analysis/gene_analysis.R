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

#read in microbe census data
census <- read_delim(file = paste(wd, "/data/microbe_census.txt", sep = ""),
                     delim = "\t", col_types = list(col_character(), 
                                                    col_number(),
                                                    col_number(), 
                                                    col_number()))

#Read in diversity data (rplB summary)
summarised <- read_delim(paste(wd, "/output/rplB.summary.scg.txt", sep = ""),
                         delim = " ", col_names = TRUE)

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_JGI_map.txt", 
                                    sep=""), sep=" ", header=TRUE))

#make colores for rarefaction curves (n=12)
rarecol <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow",
             "hotpink", "green", "red", "brown", "grey", "purple")

#prep colors for diversity
color <- c("#FF7F00", "#7570B3", "#CAB2D6", "#FBB4AE", "#F0027F", "#BEBADA", "#E78AC3", "#A6D854", "#B3B3B3", "#386CB0", "#BC80BD", "#FFFFCC", "#BF5B17", "#984EA3", "#CCCCCC", "#FFFF99", "#B15928", "#F781BF", "#FDC086", "#A6CEE3", "#FDB462", "#FED9A6", "#E6AB02", "#E31A1C", "#B2DF8A", "#377EB8", "#FCCDE5", "#80B1D3", "#FFD92F", "#33A02C", "#66C2A5", "#666666")

#make color pallette
GnYlOrRd <- colorRampPalette(colors=c("green", "yellow", "orange","red"), 
                             bias=2)

##############################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES FOR acr3#
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
  summarise(N = length(Site), Rel.Total = sum(Fraction.Abundance))

#change "cen" to "Cen" so it matches outside data
data.acr3$Site <- gsub("cen", "Cen", data.acr3$Site)

#make column for organism name and join with microbe census data and normalize to it
data.acr3 <- data.acr3 %>%
  left_join(census, by = "Site") %>%
  left_join(summarised, by = "Site") %>%
  left_join(meta, by = "Site") %>%
  separate(Taxon, into = c("coded_by", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), 
           sep = ",definition=") %>%
  rename(Temp = SoilTemperature_to10cm) %>%
  select(Site, As_ppm, Temp, Classification, organism, Abundance, 
         Fraction.Abundance, GE, rplB) %>%
  mutate(Normalized.Abundance.census = Abundance / GE, 
         Normalized.Abundance.rplB = Abundance / rplB)

#get taxonomy for organisms
acr3.ncbi <- tax_name(query = data.acr3$organism, 
                      get = c("genus", "class", "phylum"), db = "ncbi")

#change query column to "organism"
acr3.ncbi$organism <- acr3.ncbi$query

#join taxanomic information with abundance information
data.acr3.t <- data.acr3 %>%
  left_join(acr3.ncbi, by = "organism") %>%
  unique()

#get data counts
acr3.summary <- data.acr3 %>%
  mutate(Gene = "acr3", Type = "AsRG")

#save these numbers for meta analysis
write.table(acr3.summary, file = paste(wd, "/output/acr3.summary.txt", 
                                       sep = ""), row.names = FALSE)

#change NA phylum to metagenome
data.acr3.t$phylum[is.na(data.acr3.t$phylum)] = "Metagenome"


#order based on temperature
data.acr3.t$Site <- factor(data.acr3.t$Site, 
                           data.acr3.t$Site[order(data.acr3.t$Temp)])

(acr3.abundance.tax.bar.plot <- ggplot(data.acr3.t, aes(x = Site, y = Normalized.Abundance.census, fill = phylum)) +
    geom_bar(stat = "identity") +
    ylab("acr3 Abundance (normalized to genome equivalents)") +
    scale_fill_manual(values = color) +
    theme_classic(base_size = 12))
ggsave(acr3.abundance.tax.bar.plot, paste(wd, "/figures/acr3.census.phylum.eps", sep=""), width = 10)

#order based on temperature
data.acr3$Site <- factor(data.acr3$Site, 
                           levels = data.acr3$Site[order(data.acr3$Temp)])

#plot data
(acr3.abundance.plot <- ggplot(data.acr3.t, 
                               aes(x = Site, y = Normalized.Abundance.rplB,
                                   fill = phylum)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color) +
    ylab("acr3 Abundance (normalized to rplB)") +
    theme_classic())

ggsave(acr3.abundance.plot, filename = paste(wd, "/figures/acr3.rplB.phylum.png", sep=""))

(acr3.abundance.census.plot <- ggplot(data.acr3, aes(x = Site,  y = Normalized.Abundance.census, fill = Classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("firebrick2", "yellow1", "green3")) +
    ylab("acr3 Abundance (normalized to genome equivalents)") +
    theme_classic())

ggsave(acr3.abundance.census.plot, 
       filename = paste(wd, "/figures/acr3.abundance.census.png", sep=""))


(acr3.abundance.census.taxon.plot <- ggplot(data.acr3, aes(x = organism,  y = Normalized.Abundance.census)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("acr3 abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(acr3.abundance.census.taxon.plot, 
       paste(wd, "/figures/acr3.abundance.census.taxon.png", sep=""), 
       height = 10)

(acr3.abundance.taxon.plot <- ggplot(data.acr3, aes(x = organism, y = Normalized.Abundance.rplB)) +
    geom_point(aes(color = Temp, shape = Classification)) +
    ylab("acr3 abundance (normalized to rplB)") +
    xlab("Taxon") +
    coord_flip() +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(acr3.abundance.taxon.plot, paste(wd, "/figures/acr3.abundance.taxon.png", sep=""), height = 10)

#summarise acr3 information
data.acr3.sum <- data.acr3 %>%
  group_by(Classification, Site, Temp) %>%
  summarise(Abundance.census.tot = sum(Normalized.Abundance.census), 
            Abundance.rplB.tot = sum(Normalized.Abundance.rplB))

(acr3.abundance.census.box <- ggplot(data.acr3.sum, aes(x = Classification, y = Abundance.census.tot)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, aes(color = Temp)) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("acr3 abundance (normalized to genome equivalents)"))

ggsave(acr3.abundance.census.box, 
       filename = paste(wd, "/figures/acr3.abundance.census.box.png", sep=""))

(acr3.abundance.rplB.box <- ggplot(data.acr3.sum, aes(x = Classification, y = Abundance.rplB.tot)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, aes(color = Temp)) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("acr3 abundance (normalized to genome equivalents)"))

ggsave(acr3.abundance.rplB.box, 
       filename = paste(wd, "/figures/acr3.abundance.rplB.box.png", sep=""))

###########################
#acr3 DIVERSITY ANALYSIS#
###########################

#read in distance matrix
acr3 <- read.delim(file = paste(wd, "/data/acr3_rformat_dist_0.03.txt", 
                                sep=""))
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
rarecurve(otu.acr3, step=5, col = rarecol, label = TRUE)

#rarefy
rare.acr3=rarefy_even_depth(otu.acr3, sample.size = min(sample_sums(otu.acr3)), rngseed = TRUE)

#check curve
rarecurve(rare.acr3, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = TRUE)

##make biom for phyloseq
phylo.acr3=merge_phyloseq(rare.acr3, metad)

#plot phylo richness
(acr3.richness=plot_richness(phylo.acr3, x="Classification", color = "SoilTemperature_to10cm") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

#save plot 
ggsave(acr3.richness, filename = paste(wd, "/figures/acr3.richness.png", sep=""), 
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
(acr3.evenness <- ggplot(plieou.acr3, aes(x = Classification, y = plieou.acr3)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Evenness") +
    xlab("Fire classification") +
    theme_bw(base_size = 12))

#save evenness plot
ggsave(acr3.evenness, filename = paste(wd, "/figures/acr3.evenness.png", sep=""), 
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

