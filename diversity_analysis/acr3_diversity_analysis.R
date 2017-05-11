#load dependencies 
library(phyloseq)
library(ape)
library(biomformat)
library(vegan)
library(tidyverse)
library(reshape2)

#print working directory for future references
wd=print(getwd())

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_JGI_map.txt", sep=""), sep=" ", header=TRUE))

#remove Cen16 from metadata since we don't have acr3 info yet
meta=meta[which(meta$Site == "Cen14" | meta$Site == "Cen03"),]
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

##NEW
#acr3 <- data.frame(t(acr3))
#acr3 <- acr3[-1,]
#colnames(acr3) <- c("OTU", "Cen03", "Cen14")
#acr3$Cen03 <- as.numeric(as.character(acr3$Cen03))
#acr3$Cen14 <- as.numeric(as.character(acr3$Cen14))

#acr3 <- acr3 %>%
#  mutate(Cen03 = Cen03/1073, Cen14= Cen14/1839)

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
    geom_point(aes(color = SoilTemperature_to10cm), size=3) +
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
tree <- read.tree(file = paste(wd, "/data/acr3_0.03_tree2.nwk", sep=""))
tree <- phy_tree(tree)

#merge
phylo=merge_phyloseq(tree, rare, metad)

#plot tree
(tree.plot <- plot_tree(phylo, color = "SoilTemperature_to10cm", size = "abundance",
                        shape = "Classification", label.tips="taxa_names", 
                        text.size=2, ladderize="left", base.spacing = 0.03) +
    theme(legend.position = "right", legend.title = element_text(size=11),
    legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

plot_bar(phylo, fill = "Classification")

#####################################
#EXAMINE TAXON ABUNDANCE DIFFERENCES#
#####################################
#read in each file (one per site)
cen14 <- read_delim(file = paste(wd, "/data/acr3_taxonabund_cen14.txt", sep = ""), 
                    col_names = TRUE, delim = "\t")
cen03 <- read_delim(file = paste(wd, "/data/acr3_taxonabund_cen03.txt", sep = ""), 
                    col_names = TRUE, delim = "\t")

#make column for organism name
cen14 <- cen14 %>%
  mutate(Site = "Cen14", Gene = "acr3", Census = 5243, rplB = 1974.53) %>%
  mutate(Normalized.Abundance.rplB = Abundance / rplB, 
         Normalized.Abundance.census = Abundance / Census) %>%
  separate(Taxon, into = c("coded_by", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), sep = ",definition=")

cen03 <- cen03 %>%
  mutate(Site = "Cen03", Gene = "acr3", Census = 4316, rplB = 1347.908) %>%
  mutate(Normalized.Abundance.rplB = Abundance / rplB, 
         Normalized.Abundance.census = Abundance / Census) %>%  
  separate(Taxon, into = c("coded_by", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), sep = ",definition=")

#join data together
acr3 <- rbind(cen03, cen14)
acr3.high <- acr3[which(acr3$Abundance > 10),]

#plot data
ggplot(acr3, aes(x = Site, y = Normalized.Abundance.rplB)) +
  geom_bar(stat = "identity")

ggplot(acr3, aes(x = Site, y = Normalized.Abundance.census)) +
  geom_bar(stat = "identity")

ggplot(acr3, aes(x = organism, y = Normalized.Abundance.census)) +
  geom_point(aes(color = Site)) +
  coord_flip()

ggplot(acr3, aes(x = organism, y = Normalized.Abundance.rplB)) +
  geom_point(aes(color = Site)) +
  coord_flip()

  
names <- acr3$organism
write.table(names, file = "/Users/dunivint/Downloads/acr3.taxonomy.txt", row.names = FALSE, quote = FALSE)
  
  

  
  
  