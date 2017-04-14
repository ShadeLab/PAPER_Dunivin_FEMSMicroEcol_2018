#load dependencies 

#set working directory
setwd("/Users/dunivint/Documents/ShadeLab/Experiments/Xander/R_analysis/")

#read in metadata
meta.data=data.frame(read.delim(file = "meta.txt", header=TRUE))

#make sample names row names for metadata
row.names(meta.data)=meta.data[,1]
meta.data=meta.data[,-1]

#call metadata sample data
metad=sample_data(meta.data)

#read in distance matrix
rpoB=data.matrix(read.delim(file = "rformat_dist_0.03.txt", header = TRUE))

#remove first column
rpoB=rpoB[,-1]

#add row names back
rownames(rpoB)=c("cen01_rplB_45_final_prot_aligned", "cen07_rplB_45_final_prot_aligned", "cen12_rplB_45_final_prot_aligned", "cen10_rplB_45_final_prot_aligned")

#otu table
otu=otu_table(rpoB, taxa_are_rows = FALSE)

#see rarefaction curve
rarecurve(otu, step=5)

#rarefy
rare=rarefy_even_depth(otu, sample.size = min(sample_sums(otu)), rngseed = TRUE)

#check curve
rarecurve(rare, step=5)

#read tree (need to fix names)
tree=read.tree("test_rplB_tree.nwk")

##make biom for phyloseq
phylo=merge_phyloseq(rare, metad)

#plot phylo richness
plot_richness(phylo, x="Site", shape="History")

#plot ordination
ord <- ordinate(phylo, method="PCoA", distance="bray")
plot_ordination(phylo, ord, color="Site", shape="History", title="Bray Curtis") +
  geom_point(size=5) +
  theme_light(base_size = 12)

#plot tree with abundance :) 
plot_tree(phylo, color="Site", size="abundance", label.tips="OTU", text.size=3, ladderize="left") +
  theme(legend.position = "bottom", legend.title = element_text(size=12), legend.key = element_blank())




