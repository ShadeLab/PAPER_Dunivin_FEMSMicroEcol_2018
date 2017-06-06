library(vegan)
library(igraph)
library(Hmisc)
library(psych)

ClassC.t <- t(ClassC)
ClassC_norm <- ClassC.t
for(i in 1:9){ClassC_norm[,i]=ClassC.t[,i]/sum(ClassC.t[,i])}

#make presence absence matrix
ClasC_normPA <- (ClassC_norm>0)*1

#list all OTUs that are in less than half of the sites


corr <- corr.test(ClassC, method = "spearman", adjust = "fdr")
print(corr$r)
print(corr$p)
net <- graph_from_data_frame(d=corr$r, directed = T)
plot(net, edge.arrow.size=0.4, vertex.label = rownames(ClassC_norm))
