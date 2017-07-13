library(vegan)
library(psych)
library(tidyverse)
library(qgraph)
library(phyloseq)

#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis
wd <- print(getwd())

#setwd to diversity analysis
setwd(paste(wd, "/diversity_analysis", sep = ""))

#now change wd obj
wd <- print(getwd())

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_FULL_map.txt", 
                                    sep=""), sep=" ", header=TRUE))

#write OTU naming function
naming <- function(file) {
  gsub("OTU", deparse(substitute(file)), colnames(file))
}

#temporarily change working directories
setwd(paste(wd, "/data/0.1_clust", sep = ""))

#list filenames of interest
filenames <- list.files(pattern="*_rformat_dist_0.1.txt")

#move back up directories
setwd("../..")

#make dataframes of all OTU tables
for(i in filenames){
  filepath <- file.path(paste(wd, "/data/0.1_clust", sep = ""),paste(i,sep=""))
  assign(gsub("_rformat_dist_0.1.txt", "", i), read.delim(filepath,sep = "\t"))
}

#change OTU to gene name
colnames(acr3) <- naming(acr3)
colnames(aioA) <- naming(aioA)
colnames(arsB) <- naming(arsB)
colnames(`AAC6-Ia`) <- naming(`AAC6-Ia`)
colnames(adeB) <- naming(adeB)
colnames(arrA) <- naming(arrA)
colnames(arsA) <- naming(arsA)
colnames(arsC_glut) <- naming(arsC_glut)
colnames(arsC_thio) <- naming(arsC_thio)
colnames(arsD) <- naming(arsD)
colnames(arsM) <- naming(arsM)
colnames(arxA) <- naming(arxA)
colnames(CEP) <- naming(CEP)
colnames(ClassA) <- naming(ClassA)
colnames(ClassB) <- naming(ClassB)
colnames(ClassC) <- naming(ClassC)
colnames(dfra12) <- naming(dfra12)
colnames(rplB) <- naming(rplB)
colnames(intI) <- naming(intI)
colnames(sul2) <- naming(sul2)
colnames(tetA) <- naming(tetA)
colnames(tetW) <- naming(tetW)
colnames(tetX) <- naming(tetX)
colnames(tolC) <- naming(tolC)
colnames(vanA) <- naming(vanA)
colnames(vanH) <- naming(vanH)
colnames(vanX) <- naming(vanX)
colnames(vanZ) <- naming(vanZ)

#join together all files
otu_table <- acr3 %>%
  left_join(aioA, by = "X") %>%
  left_join(`AAC6-Ia`, by = "X") %>%
  left_join(arrA, by = "X") %>%
  left_join(arsA, by = "X") %>%
  left_join(arsB, by = "X") %>%
  left_join(arsC_glut, by = "X") %>%
  left_join(arsC_thio, by = "X") %>%
  left_join(arsD, by = "X") %>%
  left_join(arsM, by = "X") %>%
  left_join(arxA, by = "X") %>%
  left_join(CEP, by = "X") %>%
  left_join(ClassA, by = "X") %>%
  left_join(ClassB, by = "X") %>%
  left_join(ClassC, by = "X") %>%
  left_join(dfra12, by = "X") %>%
  left_join(intI, by = "X") %>%
  left_join(rplB, by = "X") %>%
  left_join(sul2, by = "X") %>%
  left_join(tetA, by = "X") %>%
  left_join(tetW, by = "X") %>%
  left_join(tetX, by = "X") %>%
  left_join(tolC, by = "X") %>%
  left_join(vanA, by = "X") %>%
  left_join(vanH, by = "X") %>%
  left_join(vanX, by = "X") %>%
  left_join(vanZ, by = "X") %>%
  rename(Site =X) 

#check total diversity?
#otu_table_div <- otu_table
#row.names(otu_table_div) <- otu_table_div[,1]
#otu_table_div <- otu_table_div[,-1]
#otu_table_div[is.na(otu_table_div)] <- 0
#otu <- otu_table(otu_table_div, taxa_are_rows = FALSE)
#rarecurve(otu, step=5, label = TRUE)
#plot_richness(otu, measures = c("Simpson", "Fisher"))
#a <- colSums(otu_table_div != 0)
#sum(a)
#read in rplB data
rplB <- read_delim(paste(wd, "/output/rplB.summary.scg.txt", sep = ""), delim  = " ")

#get mean rplB
mean.rplB <- mean(rplB$rplB)

#add rplB data to otu_table
otu_table.rplB <- rplB %>%
  left_join(otu_table, by = "Site") %>%
  mutate(ratio.rplB = rplB/mean.rplB)



#normalize to rplB
otu_table_norm <- otu_table.rplB
for(i in 4:5662){otu_table_norm[,i]=otu_table.rplB[,i]/otu_table.rplB[,5662]}

#add in metadata
otu_table_norm_annotated <- otu_table_norm %>%
  left_join(meta, by = "Site") %>%
  select(Site, acr3_001:ratio.rplB, As_ppm, SoilTemperature_to10cm, OrganicMatter_500:Fe_ppm) %>%
  select(-ratio.rplB)


#change to df and add row names back
otu_table_norm_annotated <- as.data.frame(otu_table_norm_annotated)
rownames(otu_table_norm_annotated) <- otu_table_norm_annotated[,1]

#remove first column
otu_table_norm_annotated=otu_table_norm_annotated[,-1]

#make data matrix
otu_table_norm_annotated=data.matrix(otu_table_norm_annotated)

#transpose data
otu_table_norm_annotated.t <- t(otu_table_norm_annotated)

#replace NAs with zeros
otu_table_norm_annotated.t[is.na(otu_table_norm_annotated.t)] <- 0


#make presence absence matrix
otu_table_normPA <- (otu_table_norm_annotated.t>0)*1

#list OTUs present in less than 2 samples
abund <- otu_table_normPA[which(rowSums(otu_table_normPA) > 4),]

#remove OTUs with presence in less than 4 sites
otu_table_norm.slim <- otu_table_norm_annotated.t[which(rownames(otu_table_norm_annotated.t) %in%
                              rownames(abund)),]

otu_table_norm.export <- otu_table_norm_annotated.t
#replace all 0's with NA for export
is.na(otu_table_norm.export) <- !otu_table_norm.export

#replace NAs with nothing
otu_table_norm.export[is.na(otu_table_norm.export)] <- ""

#write table as output for SparCC
write.table(otu_table_norm.export, 
            file = (paste(wd, "/output/otu_table.txt", 
                          sep = "")), 
            sep = "\t", quote = FALSE)

#transpose dataset
otu_table_norm.slim.t <- t(otu_table_norm.slim)


#find correlations between OTUs
corr <- corr.test(otu_table_norm.slim.t, 
                  method = "spearman", adjust = "fdr", alpha = 0.01)

corr.r <- as.matrix(print(corr$r, long = TRUE))
corr.p <- as.matrix(print(corr$p, long = TRUE))

corr.r[which(corr.p > 0.01)] <- 0

#save correlation data
write.table(corr$r, paste(wd, "/output/corr_table.rplB.txt", sep = ""), quote = FALSE)

#make network of correlations
qgraph(corr$r, minimum = "sig", sampleSize=13, 
       layout = "spring", details = TRUE,
       graph = "cor", label.cex = 3, curve = 0.2, curveAll = TRUE,
       alpha = 0.01)


library(igraph)
cor_mat<-as.matrix(corr.r)
diag(cor_mat)<-0
graph<-graph.adjacency(cor_mat,weighted=TRUE,mode="lower")
graph <- delete.edges(graph, E(graph)[ abs(weight) < 0.68])
#graph <- delete.vertices(graph,which(degree(graph)<1))
E(graph)$color <- "grey";
E(graph)$width <- 1;
E(graph)[weight > 0]$color <- "green";
E(graph)[weight < 0]$color <- "red";
V(graph)$color <- "grey";
E(graph)$width <- abs(E(graph)$weight*2);
tkplot(graph, edge.curved = FALSE)

plot(graph)
get.edgelist(graph)

#test without rplB!
#remove column based on pattern (rplB)
otu_table_norm.slim.t.genes <- otu_table_norm.slim.t[, -grep("rplB", colnames(otu_table_norm.slim.t))]

corr.genes <- corr.test(otu_table_norm.slim.t.genes, 
                  method = "spearman", adjust = "fdr", alpha = 0.01)

corr.r.genes <- as.matrix(print(corr.genes$r, long = TRUE))
corr.p.genes <- as.matrix(print(corr.genes$p, long = TRUE))

corr.r.genes[which(corr.p.genes > 0.01)] <- 0

library(igraph)
cor_mat<-as.matrix(corr.r.genes)
diag(cor_mat)<-0
graph<-graph.adjacency(cor_mat,weighted=TRUE,mode="lower")
graph <- delete.edges(graph, E(graph)[ abs(weight) < 0.68])
#graph <- delete.vertices(graph,which(degree(graph)<1))
E(graph)$color <- "grey";
E(graph)$width <- 1;
E(graph)$width <- E(graph)$weight*2;
E(graph)[weight > 0]$color <- "green";
E(graph)[weight < 0]$color <- "red";
V(graph)$color <- "grey";
tkplot(graph, edge.curved = FALSE)
