#load dependencies 
library(phyloseq)
library(vegan)
library(ggplot2)

#read in metadata
meta=data.frame(read.delim(file = "Centralia_full_map.txt", sep=" ", header=TRUE))

#read in distance matrix
rplB=read.delim(file = "rformat_dist_0.03.txt")

#make sample names row names for metadata
row.names(meta)=meta[,1]
meta=meta[,-1]

#call metadata sample data
metad=sample_data(meta)

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
rarecurve(rare, step=5, col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), label = FALSE)

##make biom for phyloseq
phylo=merge_phyloseq(rare, metad)

#plot phylo richness
plot_richness(phylo, x="Site", shape="History", color = "Temperature")

#plot ordination
ord <- ordinate(phylo, method="PCoA", distance="bray")
plot_ordination(phylo, ord, color="Site", shape="History", title="Bray Curtis") +
  geom_point(size=5) +
  theme_light(base_size = 12)

#plot nmds
ord1 <- ordinate(phylo, method="NMDS", distance="bray")
plot_ordination(phylo, ord1, color="Site", shape="History") +
  geom_point(size=5)

