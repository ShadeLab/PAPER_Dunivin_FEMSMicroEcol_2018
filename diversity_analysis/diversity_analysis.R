#load dependencies 
library(phyloseq)
library(vegan)
library(tidyverse)

#print working directory for future references
wd=print(getwd())

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
  geom_jitter(aes(color = SoilTemperature_to10cm), size=5) +
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
  theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord, filename = paste(wd, "/figures/bc.ord.png", sep=""), 
       width = 6, height = 5)


#plot Sorenson ordination
ord.sor <- ordinate(phylo, method="PCoA", distance="bray", binary = TRUE)
(sorenson.ord=plot_ordination(phylo, ord.sor, shape="Classification", title="Sorenson") +
    geom_jitter(aes(color = SoilTemperature_to10cm), size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save bray curtis ordination
ggsave(sorenson.ord, filename = paste(wd, "/figures/bc.ord.png", sep=""), 
       width = 6, height = 5)


#make object phylo with tree and biom info
tree <- read.tree(file = paste(wd, "/data/rplB_0.03_tree.nwk", sep=""))
tree <- phy_tree(tree)

#merge
phylo=merge_phyloseq(tree, rare, metad)

#plot tree
(tree.plot <- plot_tree(phylo, shape = "Classification", size = "abundance",
                        color = "SoilTemperature_to10cm", label.tips=NULL, 
                        text.size=2, ladderize="left", base.spacing = 0.04) +
    theme(legend.position = "right", legend.title = element_text(size=11),
          legend.key =element_blank()) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")))

ggsave(tree.plot, filename = paste(wd, "/figures/rplb.tree.png", sep=""), 
       width = 10, height = 40)

