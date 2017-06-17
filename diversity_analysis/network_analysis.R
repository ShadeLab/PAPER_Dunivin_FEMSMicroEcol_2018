library(vegan)
library(psych)
library(tidyverse)
library(qgraph)

#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis
wd <- print(getwd())

#setwd to diversity analysis
setwd(paste(wd, "/diversity_analysis", sep = ""))

#now change wd obj
wd <- print(getwd())

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_JGI_map.txt", 
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
colnames(acr3) <- gsub("OTU", deparse(substitute(acr3)), colnames(acr3))
colnames(aioA) <- gsub("OTU", deparse(substitute(aioA)), colnames(aioA))
colnames(arsC_glut) <- gsub("OTU", deparse(substitute(arsC_glut)), colnames(arsC_glut))
colnames(arsC_thio) <- gsub("OTU", deparse(substitute(arsC_thio)), colnames(arsC_thio))
colnames(arsM) <- gsub("OTU", deparse(substitute(arsM)), colnames(arsM))
colnames(ClassA) <- gsub("OTU", deparse(substitute(ClassA)), colnames(ClassA))
colnames(ClassB) <- gsub("OTU", deparse(substitute(ClassB)), colnames(ClassB))
colnames(ClassC) <- gsub("OTU", deparse(substitute(ClassC)), colnames(ClassC))
colnames(intI) <- gsub("OTU", deparse(substitute(intI)), colnames(intI))
colnames(vanA) <- gsub("OTU", deparse(substitute(vanA)), colnames(vanA))
colnames(vanX) <- gsub("OTU", deparse(substitute(vanX)), colnames(vanX))
colnames(vanZ) <- gsub("OTU", deparse(substitute(vanZ)), colnames(vanZ))
colnames(vanH) <- gsub("OTU", deparse(substitute(vanH)), colnames(vanH))
colnames(rplB) <- gsub("OTU", deparse(substitute(rplB)), colnames(rplB))

#join together all files
otu_table <- acr3 %>%
  left_join(aioA, by = "X") %>%
  left_join(ClassA, by = "X") %>%
  left_join(ClassB, by = "X") %>%
  left_join(ClassC, by = "X") %>%
  left_join(intI, by = "X") %>%
  left_join(vanA, by = "X") %>%
  left_join(vanX, by = "X") %>%
  left_join(vanZ, by = "X") %>%
  left_join(vanH, by = "X") %>%
  left_join(arsC_glut, by = "X") %>%
  left_join(arsC_thio, by = "X") %>%
  left_join(rplB, by = "X") %>%
  rename(Site = X)

#count otu columns
otus <- ncol(otu_table)

#add in metadata
otu_table <- otu_table %>%
  left_join(meta, by = "Site") %>%
  select(Site:As_ppm, SoilTemperature_to10cm, OrganicMatter_500:Fe_ppm)

#add row names back
rownames(otu_table)=otu_table[,1]

#remove first column
otu_table=otu_table[,-1]

#make data matrix
otu_table=data.matrix(otu_table)

#replace NAs with zeros
otu_table[is.na(otu_table)] <- 0

otu_table.t <- t(otu_table)
otu_table_norm <- otu_table.t

for(i in 1:12){otu_table_norm[,i]=otu_table.t[,i]/sum(otu_table.t[,i])}

#make presence absence matrix
otu_table_normPA <- (otu_table_norm>0)*1

#list OTUs present in less than half of samples
abund <- otu_table_normPA[which(rowSums(otu_table_normPA) > 3),]

#remove OTUs with presence in less than 4 sites
otu_table_norm.slim <- otu_table_norm[which(rownames(otu_table_norm) %in%
                              rownames(abund)),]

#transpose dataset
otu_table_norm.slim.t <- t(otu_table_norm.slim)

#find correlations between OTUs
corr <- corr.test(otu_table_norm.slim.t, method = "spearman", adjust = "fdr")

#make network of correlations
qgraph(corr$r, minimum = "sig", sampleSize=12, layout = "spring", details = TRUE, graph = "cor", label.cex = 2, alpha = 0.01)
