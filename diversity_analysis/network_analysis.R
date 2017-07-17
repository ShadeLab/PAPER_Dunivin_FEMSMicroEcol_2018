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
  left_join(adeB, by = "X") %>%
  rename(Site =X) 

#list otus that are gene matches
mixed <- c("aioA_13", "aioA_31", "aioA_34", "aioA_36", "aioA_37", "aioA_42", "aioA_43", "aioA_49", "aioA_51", "arrA_1", "arrA_2", "arrA_3", "arrA_4", "arrA_5", "arrA_6", "arrA_8","arxA_11", "arxA_04")

#remove OTUs that are gene matches
otu_table <- otu_table[,!names(otu_table) %in% mixed]

#rename arxA column that is actually arrA (with no match)
otu_table <- otu_table %>%
  rename(arrA_3 = arxA_03)

#write file to save OTU table
write.table(otu_table, paste(wd, "/output/otu_table.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

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
for(i in 4:6572){otu_table_norm[,i]=otu_table.rplB[,i]/otu_table.rplB[,3]}

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



############################
#MAKE GENE ABUNDANCE GRAPHS#
############################
#read in gene classification data
gene <- read_delim(paste(wd, "/data/gene_classification.txt",  sep=""), 
                   delim = "\t", col_names = TRUE)

#melt data to separate otu and abundance per site
gene_abundance <- otu_table_norm %>%
  melt(variable.names = c("Site", "OTU"), value.name = "RelativeAbundance") %>%
  rename(OTU = variable)

#change arsC_ to arsC
gene_abundance$OTU <- gsub("arsC_", "arsC", gene_abundance$OTU)

gene_abundance <- gene_abundance %>%
  separate(OTU, into = "Gene", sep = "_", remove = FALSE) %>%
  left_join(gene, by = "Gene") %>%
  left_join(meta, by = "Site") %>%
  unique()

#replace NAs with zeros
gene_abundance$RelativeAbundance[is.na(gene_abundance$RelativeAbundance)] = 0

#remove total rows
gene_abundance <- gene_abundance[-which(gene_abundance$OTU == "Total"),]
gene_abundance <- gene_abundance[-which(gene_abundance$OTU == "rplB"),]
gene_abundance <- gene_abundance[-which(gene_abundance$OTU == "ratio.rplB"),]

#make color pallette for Centralia temperatures
GnYlOrRd <- colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)


#prep colors for phylum diversity
color <- c("#FF7F00", "#7570B3", "#CAB2D6", "#FBB4AE", "#F0027F", "#BEBADA", "#E78AC3", "#A6D854", "#B3B3B3", "#386CB0", "#BC80BD", "#FFFFCC", "#BF5B17", "#984EA3", "#CCCCCC", "#FFFF99", "#B15928", "#F781BF", "#FDC086", "#A6CEE3", "#FDB462", "#FED9A6", "#E6AB02", "#E31A1C", "#B2DF8A", "#377EB8", "#FCCDE5", "#80B1D3", "#FFD92F", "#33A02C", "#66C2A5", "#666666", "black", "brown", "grey", "red", "blue", "green", "orange")


#summarise data
gene_abundance_summary <- gene_abundance %>%
  group_by(Group, Description, Gene, Classification, Site, SoilTemperature_to10cm) %>%
  summarise(Total = sum(RelativeAbundance))

#plot!
(barplot <- ggplot(gene_abundance_summary, aes(x = Classification, 
                                                  y = Total)) +
    geom_bar(stat = "identity", aes(fill = Gene)) +
    facet_wrap(~Group) +
    scale_fill_manual(values = color) +
    ylab("Total gene count (normalized to rplB)") +
    theme_classic(base_size = 15))

(boxplot <- ggplot(gene_abundance_summary, aes(x = Classification, 
                                               y = Total)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm)) +
    facet_wrap(~Gene, scales = "free_y") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                         guide_legend(title="Temperature (°C)")) +
    ylab("Total gene count (normalized to rplB)") +
    theme_classic(base_size = 12))

######################
#CORRELATION ANALYSES#
######################

#decast for abundance check but first make Site a character
gene_abundance_summary$Site <- as.character(gene_abundance_summary$Site)
cast.gene <- data.frame(acast(gene_abundance_summary, Site ~ Gene, value.var = "Total"))

#call na's zeros
cast.gene[is.na(cast.gene)] =0

#remove unnecessary sites in meta data
meta <- meta[which(meta$Site %in% otu_table$Site),]

#correlate data with temp
corr <- print(corr.test(x = cast.gene, y = meta[,c(2,16, 20:30)], method = "spearman", adjust = "fdr"), short = FALSE)
cor <- corr.test(x = cast.gene, y = meta[,c(2,16, 20:30)], method = "spearman", adjust = "fdr")
#arsenic and antibiotic resistance genes are not correlated with temperature!

#check correlations between genes
corr.genes <- print(corr.test(cast.gene, method = "spearman", adjust = "fdr"), 
                    short = FALSE)
corr.genes.matrix <- corr.test(cast.gene, method = "spearman", adjust = "fdr")
cor.plot(t(corr.genes.matrix$r), stars = TRUE, pval=corr.genes.matrix$p, numbers = TRUE, diag = FALSE, xlas=2)


######################
#MANN WHITNEY U TESTS#
######################

#remove C17 (testing recovered v reference)
data.mann <- gene_abundance_summary[-which(gene_abundance_summary$Classification == "Reference"),]

#subset data for each gene to test:
#acr3
acr3 <- subset(x = data.mann, subset = Gene == "acr3")
acr3.cast <- print(wilcox.test(acr3$Total~acr3$Classification, paired = FALSE))
t.test(acr3$Total ~ acr3$Classification)


#aioA
aioA <- subset(x = data.mann, subset = Gene == "aioA")
aioA.cast <- print(wilcox.test(aioA$Total~aioA$Classification, paired = FALSE))
t.test(aioA$Total ~ aioA$Classification)


#arsM
arsM <- subset(x = data.mann, subset = Gene == "arsM")
arsM.cast <- print(wilcox.test(arsM$Total~arsM$Classification, paired = FALSE))
t.test(arsM$Total ~ arsM$Classification)


#arsCglut
arsCglut <- subset(x = data.mann, subset = Gene == "arsCglut")
arsCglut.cast <- print(wilcox.test(arsCglut$Total~arsCglut$Classification, paired = FALSE))
t.test(arsCglut$Total ~ arsCglut$Classification)


#arsCthio
arsCthio <- subset(x = data.mann, subset = Gene == "arsCthio")
arsCthio.cast <- print(wilcox.test(arsCthio$Total~arsCthio$Classification, paired = FALSE))
t.test(arsCthio$Total ~ arsCthio$Classification)

#arsA
arsA <- subset(x = data.mann, subset = Gene == "arsA")
arsA.cast <- print(wilcox.test(arsA$Total~arsA$Classification, paired = FALSE))
t.test(arsA$Total ~ arsA$Classification)

#arsD
arsD <- subset(x = data.mann, subset = Gene == "arsD")
arsD.cast <- print(wilcox.test(arsD$Total~arsD$Classification, paired = FALSE))
t.test(arsD$Total ~ arsD$Classification)

#intI
intI <- subset(x = data.mann, subset = Gene == "intI")
intI.cast <- print(wilcox.test(intI$Total~intI$Classification, paired = FALSE))
t.test(intI$Total ~ intI$Classification)

#ClassA
ClassA <- subset(x = data.mann, subset = Gene == "ClassA")
ClassA.cast <- print(wilcox.test(ClassA$Total~ClassA$Classification, paired = FALSE))
t.test(ClassA$Total ~ ClassA$Classification)

#ClassB
ClassB <- subset(x = data.mann, subset = Gene == "ClassB")
ClassB.cast <- print(wilcox.test(ClassB$Total~ClassB$Classification, paired = FALSE))
t.test(ClassB$Total ~ ClassB$Classification)

#ClassC
ClassC <- subset(x = data.mann, subset = Gene == "ClassC")
ClassC.cast <- print(wilcox.test(ClassC$Total~ClassC$Classification, paired = FALSE))
t.test(ClassC$Total ~ ClassC$Classification)

#vanA
vanA <- subset(x = data.mann, subset = Gene == "vanA")
vanA.cast <- print(wilcox.test(vanA$Total~vanA$Classification, paired = FALSE))
t.test(vanA$Total ~ vanA$Classification)

#vanH
vanH <- subset(x = data.mann, subset = Gene == "vanH")
vanH.cast <- print(wilcox.test(vanH$Total~vanH$Classification, paired = FALSE))
t.test(vanH$Total ~ vanH$Classification)

#vanX
vanX <- subset(x = data.mann, subset = Gene == "vanX")
vanX.cast <- print(wilcox.test(vanX$Total~vanX$Classification, paired = FALSE))
t.test(vanX$Total ~ vanX$Classification)

#vanZ
vanZ <- subset(x = data.mann, subset = Gene == "vanZ")
vanZ.cast <- print(wilcox.test(vanZ$Total~vanZ$Classification, paired = FALSE))
t.test(vanZ$Total ~ vanZ$Classification)

#tolC
tolC <- subset(x = data.mann, subset = Gene == "tolC")
tolC.cast <- print(wilcox.test(tolC$Total~tolC$Classification, paired = FALSE))
t.test(tolC$Total ~ tolC$Classification)

#sul2
sul2 <- subset(x = data.mann, subset = Gene == "sul2")
sul2.cast <- print(wilcox.test(sul2$Total~sul2$Classification, paired = FALSE))
t.test(sul2$Total ~ sul2$Classification)

#dfra12
dfra12 <- subset(x = data.mann, subset = Gene == "dfra12")
dfra12.cast <- print(wilcox.test(dfra12$Total~dfra12$Classification, paired = FALSE))
t.test(dfra12$Total ~ dfra12$Classification)

#adeB
adeB <- subset(x = data.mann, subset = Gene == "adeB")
adeB.cast <- print(wilcox.test(adeB$Total~adeB$Classification, paired = FALSE))
t.test(adeB$Total ~ adeB$Classification)


