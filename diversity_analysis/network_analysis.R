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

#read in rplB data
rplB <- read_delim(paste(wd, "/output/rplB.summary.scg.txt", sep = ""), delim  = " ")

#add rplB data to otu_table
otu_table.rplB <- rplB %>%
  left_join(otu_table, by = "Site")

#count otu columns
otus <- ncol(otu_table) -1

#normalize to rplB
otu_table_norm <- otu_table.rplB
for(i in 4:5481){otu_table_norm[,i]=otu_table.rplB[,i]/otu_table.rplB[,3]}

#add in metadata
otu_table_norm_annotated <- otu_table_norm %>%
  left_join(meta, by = "Site") %>%
  select(Site, acr3_001:As_ppm, SoilTemperature_to10cm, OrganicMatter_500:Fe_ppm)

#change to df and add row names back
otu_table_norm_annotated <- as.data.frame(otu_table_norm_annotated)
rownames(otu_table_norm_annotated) <- otu_table_norm_annotated[,1]

#remove first column
otu_table_norm_annotated=otu_table_norm_annotated[,-1]

#make data matrix
otu_table_norm_annotated=data.matrix(otu_table_norm_annotated)

#replace NAs with zeros
otu_table_norm_annotated[is.na(otu_table_norm_annotated)] <- 0

otu_table_norm_annotated.t <- t(otu_table_norm_annotated)

#write table as output for SparCC
write.table(otu_table_norm_annotated.t, 
            file = (paste(wd, "/output/otu_table_rplBn.txt", 
                          sep = "")), 
            sep = "\t", quote = FALSE)

#make presence absence matrix
otu_table_normPA <- (otu_table_norm_annotated.t>0)*1

#list OTUs present in less than half of samples
abund <- otu_table_normPA[which(rowSums(otu_table_normPA) > 7),]

#remove OTUs with presence in less than 4 sites
otu_table_norm.slim <- otu_table_norm_annotated.t[which(rownames(otu_table_norm_annotated.t) %in%
                              rownames(abund)),]

#transpose dataset
otu_table_norm.slim.t <- t(otu_table_norm.slim)

#save otu data
write.table(otu_table_norm.slim.t, paste(wd, "/output/otu_table.rplB.txt", sep = ""), quote = FALSE, row.names = TRUE)

#find correlations between OTUs
corr <- corr.test(otu_table_norm.slim.t, 
                  method = "spearman", adjust = "fdr")

#save correlation data
write.table(corr$r, paste(wd, "/output/corr_table.rplB.txt", sep = ""), quote = FALSE)

#make network of correlations
x <- qgraph(corr$r, minimum = "sig", sampleSize=12, 
       layout = "spring", details = TRUE,
       graph = "cor", label.cex = 2, 
       threshold = "fdr", curve = 1, curveAll = TRUE,
       alpha = 0.01)
centrality_auto(x)
EBICglasso(corr$r, 13)
FDRnetwork(corr$r)
findGraph(corr$r, 13)
centralityPlot(corr$r)
ggmFit(corr$r)
