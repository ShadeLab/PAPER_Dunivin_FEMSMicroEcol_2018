#load dependencies 
library(phyloseq)
library(vegan)
library(ggplot2)

#print working directory for future references
wd=print(getwd())

#read in metadata
meta=data.frame(read.delim(file = paste(wd, "/data/Centralia_mini_map.txt", sep=""), sep=" ", header=TRUE))

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
ggsave(richness, filename = paste(wd, "/figures/richness.eps", sep=""), width = 10, height = 3)

#calculate evenness
s=specnumber(rare)
h=diversity(rare, index="shannon")
plieou=h/log(s)

#save evenness number
write.table(plieou, file = paste(wd, "/output/evenness.txt", sep=""))

#plot ordination
ord <- ordinate(phylo, method="PCoA", distance="bray")
(bc.ord=plot_ordination(phylo, ord, color="Sample", shape="Classification", title="Bray Curtis") +
  geom_jitter(size=5) +
  theme_light(base_size = 12))

#save bray curtis ordination
ggsave(bc.ord, filename = paste(wd, "/figures/bc.ord.eps", sep=""), width = 5, height = 4)


#plot nmds (so far too few points)
ord1 <- ordinate(phylo, method="NMDS", distance="bray")
plot_ordination(phylo, ord1, color="Sample", shape="Classification") +
  geom_point(size=5)

#make an output of total gene count per site
gcounts=rowSums(rplB)

#########################
Phylogeny & Abundance
#########################




