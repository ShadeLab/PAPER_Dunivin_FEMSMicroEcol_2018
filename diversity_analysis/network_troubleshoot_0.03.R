library(vegan)
library(psych)
library(tidyverse)
library(qgraph)
library(phyloseq)
library(reshape2)
library(broom)
library(taxize)
library(stringr)

#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis
wd <- print(getwd())

#setwd to diversity analysis
setwd(paste(wd, "/diversity_analysis", sep = ""))
wd <- print(getwd())

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_FULL_map.txt", 
                                    sep=""), sep=" ", header=TRUE))

#####################
#SET UP CONTIG-TABLE#
#####################

#write OTU naming function
naming <- function(file) {
  gsub("OTU", deparse(substitute(file)), colnames(file))
}

#temporarily change working directories
setwd(paste(wd, "/data", sep = ""))

#list filenames of interest
filenames <- list.files(pattern="*_rformat_dist_0.03.txt")

#move back up directories
setwd("../..")

#make dataframes of all OTU tables
for(i in filenames){
  filepath <- file.path(paste(wd, "/data", sep = ""),paste(i,sep=""))
  assign(gsub("_rformat_dist_0.03.txt", "", i), read.delim(filepath,sep = "\t"))
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
  left_join(aioA, by = "X") %>% left_join(`AAC6-Ia`, by = "X") %>% left_join(arrA, by = "X") %>% left_join(arsA, by = "X") %>% left_join(arsB, by = "X") %>% left_join(arsC_glut, by = "X") %>% left_join(arsC_thio, by = "X") %>% left_join(arsD, by = "X") %>% left_join(arsM, by = "X") %>% left_join(arxA, by = "X") %>% left_join(CEP, by = "X") %>% left_join(ClassA, by = "X") %>% left_join(ClassB, by = "X") %>% left_join(ClassC, by = "X") %>% left_join(dfra12, by = "X") %>% left_join(intI, by = "X") %>% left_join(rplB, by = "X") %>% left_join(sul2, by = "X") %>% left_join(tetA, by = "X") %>% left_join(tetW, by = "X") %>% left_join(tetX, by = "X") %>% left_join(tolC, by = "X") %>% left_join(vanA, by = "X") %>% left_join(vanH, by = "X") %>% left_join(vanX, by = "X") %>% left_join(vanZ, by = "X") %>% left_join(adeB, by = "X") %>% rename(Site =X) 

#replace all NAs (from join) with zeros
otu_table[is.na(otu_table)] <- 0

#write file to save OTU table
write.table(otu_table, paste(wd, "/output/otu_table_0.03.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

#########################################
#EXTRACT AND TIDY RPLB FOR NORMALIZATION#
#########################################

#get rplB otu data
rplB <- otu_table[,grepl("rplB", names(otu_table))]

#add site names to rplB
rownames(rplB) <-  otu_table$Site

#summarise rplB by adding all OTU counts 
# in each row (total rplB/site)
rplB_summary <- data.frame(rowSums(rplB))

#make site a column in rplB
rplB_summary$Site <- rownames(rplB_summary) 

#save rplB sums
rplB_summary_save <- rplB_summary
rplB_summary_save$Site <- gsub("cen", "Cen", rplB_summary_save$Site)
write.table(rplB_summary_save, file = paste(wd, "/output/rplB.summary.scg_0.03.txt", sep = ""), row.names = FALSE)

###################################
#CONTIG-TABLE NORMALIZATION (RPLB)#
###################################

#add rplB data to otu_table
otu_table.rplB <- rplB_summary %>%
  left_join(otu_table, by = "Site") %>%
  rename(rplB = rowSums.rplB.)

#normalize to rplB
otu_table_norm <- otu_table.rplB
for(i in 3:ncol(otu_table_norm)){otu_table_norm[,i]=otu_table.rplB[,i]/otu_table.rplB[,1]}

#add in metadata
otu_table_norm$Site <- gsub("cen", "Cen", otu_table_norm$Site)
otu_table_norm_annotated <- otu_table_norm %>%
  left_join(meta, by = "Site") %>%
  mutate(DateSince_Fire = 2014-DateFire_Elick2011) %>%
  select(rplB:As_ppm, SoilTemperature_to10cm, DateSince_Fire,
         OrganicMatter_500:Fe_ppm)

#change to df and add row names back
#remove first two columns (redundant)
otu_table_norm_annotated <- as.data.frame(otu_table_norm_annotated)
rownames(otu_table_norm_annotated) <- otu_table_norm_annotated[,2]
otu_table_norm_annotated=otu_table_norm_annotated[,-c(1,2)]

#make data matrix and transpose 
otu_table_norm_annotated.t <- t(data.matrix(otu_table_norm_annotated))

#replace NAs with zeros (except date since fire)
n <- nrow(otu_table_norm_annotated.t)-11
otu_table_norm_annotated.t[1:n,][is.na(otu_table_norm_annotated.t[1:n,])] <- 0

#make presence absence matrix
otu_table_normPA <- (otu_table_norm_annotated.t>0)*1

#replace NA with 0 (date since fire)
otu_table_normPA[is.na(otu_table_normPA)] <- 0

#########################################
#RESISTANCE GENE (FULL) NETWORK ANALYSIS#
#########################################

#determine sample cutoff
cutoff <- data.frame(Number_of_Sites = rowSums(otu_table_normPA))
ggplot(cutoff, aes(x = Number_of_Sites)) +
  geom_histogram(binwidth = 1)
quantile(cutoff$Number_of_Sites, 0.96)

#list OTUs present in less than 2 samples
abund <- otu_table_normPA[which(rowSums(otu_table_normPA) > 2),]

#remove OTUs with presence in less than 4 sites
otu_table_norm.slim.t <- t(otu_table_norm_annotated.t[which(rownames(otu_table_norm_annotated.t) %in% rownames(abund)),])

#find correlations between contigs!
corr.genes <- corr.test(otu_table_norm.slim.t, method = "spearman", adjust = "fdr", alpha = 0.01)

## prepare network graphics
#read in gene classification data
gene <- read_delim(paste(wd, "/data/gene_classification.txt",  sep=""), 
                   delim = "\t", col_names = TRUE)

#Create list of shapes and colors for network
r <- data.frame(corr.genes$r)
r$gene <- rownames(r)
r$gene <- gsub("arsC_", "arsC", r$gene)
r <- r %>%
  separate(gene, c("Gene", "Number"), by = "_", remove = FALSE) %>%
  left_join(gene, by = "Gene") %>%
  select(Gene, Number, Group, gene, gene.color) %>%
  rename(OTU = gene) %>%
  mutate(Shape = "circle")

#label r appropriately
r$Group[r$Gene == "rplB"] <- "Organism"
r$Shape[r$Group == "Organism"] <- "square"
r$Group[is.na(r$Group)] <- "Metadata"
r$Shape[r$Group == "Metadata"] <- "diamond"
r$Shape[r$Group == "AntibioticResistance"] <- "triangle"

#make vector of shapes for network
shapes <- as.vector(r$Shape)

#make rplB and metadata gene colors
r$gene.color[r$Gene == "rplB"] <- "#fffac8"
r$gene.color[r$Group == "Metadata"] <- "#D3D3D3"

#read in BLAST output results
setwd(paste(wd, "/../networks/data", sep = ""))
blast.names <- list.files(pattern="*0.03.txt")
blast <- do.call(rbind, lapply(blast.names, function(X) {
  data.frame(id = basename(X), read_delim(X, delim = "\t", col_names = FALSE))}))
setwd(wd)
#tidy blast results
blast$id <- gsub("arsC_", "arsC", blast$id)
blast.tidy <- blast %>%
  separate(id, into = c("results", "Gene", "Clust"), sep = "_") %>%
  separate(X3, into = c("accno", "other"), sep = " coded_by=") %>%
  separate(other, into = c("list", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), sep = ",definition=") %>%
  rename(OTU = X1, evalue = X4, percid = X5) %>%
  select(-c(results, X2))

#remove all blast.tidy results that have e.values > 10-5
blast.tidy <- blast.tidy[blast.tidy$evalue < 0.00001,]

#subset rplB results (pre-classified)
blast.rplB <- subset(blast.tidy, Gene == "rplB")
blast.all <- subset(blast.tidy, Gene != "rplB")

#remove all rows with poor taxanomic identifiers (ie unknown or metagenome)
blast.all <- blast.all[!grepl("metagenome", blast.all$organism),]
blast.all <- blast.all[!grepl("uncultured .", blast.all$organism),]

#order blast.all by e value (lowest first)
blast.all <- blast.all[order(blast.all$evalue, decreasing = FALSE),]

#read in taxize information
taxize <- read.delim(paste(wd, "/output/taxize_0.1_results.txt", sep = ""), sep = " ")

#remove duplicate matches
blast.all.unique <- blast.all[!duplicated(blast.all[,c(1,3)]),]
blast.all.unique <- blast.all.unique %>%
  select(Gene, OTU, organism) %>%
  rename(query = organism) %>%
  left_join(taxize, by = "query")

blast.all.unique.complete <- blast.all.unique[!is.na(blast.all.unique$phylum),]  
blast.all.unique.na <- blast.all.unique[is.na(blast.all.unique$phylum),-c(4:7)] 

#fill in blast results
blast.all.unique.na[,4:8] <- tax_name(blast.all.unique.na$query, get = c("genus", "class", "phylum"), db = "ncbi")

#remove duplicate query column and merge files
blast.all.unique.na <- blast.all.unique.na[,-5]
blast.all.unique.final <- rbind(blast.all.unique.complete, blast.all.unique.na)

#save taxize results to output so that ncbi does not need to be re-searched
write.table(blast.all.unique.final, paste(wd, "/output/taxize_0.03_results.txt", sep = ""), col.names = TRUE, row.names = FALSE)

#select columns of interest in blast.all
blast.all.unique.final <- blast.all.unique.final %>%
  select(Gene, OTU, phylum, class, genus)

#tidy rplB information (ie get phyla)
blast.rplB.tidy <- blast.rplB %>%
  separate(accno, into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";") %>%
  select(Gene, OTU, phylum, class, genus, evalue)

#remove rows that have unknown phylum
blast.rplB.tidy <- blast.rplB.tidy[!grepl("metagenomes", blast.rplB.tidy$phylum),]
blast.rplB.tidy <- blast.rplB.tidy[!grepl("environmentalsamples", blast.rplB.tidy$phylum),]
blast.rplB.tidy <- blast.rplB.tidy[!grepl("artificialsequences", blast.rplB.tidy$phylum),]

#order blast.all by e value (lowest first) and then remove duplicate rows
blast.rplB.tidy <- blast.rplB.tidy[order(blast.rplB.tidy$evalue, decreasing = FALSE),]
blast.rplB.tidy.unique <- blast.rplB.tidy[!duplicated(blast.rplB.tidy$OTU) & duplicated(blast.rplB.tidy$Gene),]

#select necessary columns in rplB tidy
blast.rplB.tidy.unique <- select(blast.rplB.tidy.unique, Gene, OTU, phylum, class, genus)

#join together rplB and blast data 
blast.final <- rbind(blast.rplB.tidy.unique, blast.all.unique.final)

#replace "OTU" with gene name so it can be joined with shape
blast.final$OTU <- with(blast.final, str_replace_all(blast.final$OTU, fixed("OTU"), blast.final$Gene))

#save taxize results to output so that ncbi does not need to be re-searched
write.table(blast.final, paste(wd, "/output/taxize_0.03_results_FINAL.txt", sep = ""), col.names = TRUE, row.names = FALSE)

#add leading zeros to r otus
r.final <- r %>%
  separate(OTU, into = c("beginning", "number"), sep = "_")

r.final$number <- sprintf("%04s", r.final$number)
r.final$OTU <- paste(r.final$beginning, r.final$number, sep = "_")

#join aesthetic (shape, taxonomy) information
aesthetics <- r.final %>%
  left_join(blast.final, by = c("OTU", "Gene")) %>%
  unique()

#read in colors for phyla
phylum.colors <- read_delim(paste(wd, "/../networks/data/phylum_colors.txt", sep = ""), delim = " ", col_names = c("phy.color", "phylum"))

#add phylum colors to aesthetics
aesthetics <- aesthetics %>%
  left_join(phylum.colors, by = "phylum")

#make taxanomic groups for network
phylum.group <- split(as.numeric(rownames(aesthetics)), list(aesthetics$phylum)) 
class.group <- split(as.numeric(rownames(aesthetics)), list(aesthetics$class), drop=TRUE) 
genus.group <- split(as.numeric(rownames(aesthetics)), list(aesthetics$genus), drop=TRUE) 

#make network of correlations
clust.network <- qgraph(corr.genes$r, minimum = "sig", sampleSize=13, 
                        details = TRUE, layout = "spring",
                        graph = "cor",label.cex = 0.5,
                        alpha = 0.01, graph = "fdr", labels = aesthetics$Gene,  label.scale.equal = TRUE, label.scale = FALSE,shape = aesthetics$Shape, node.resolution = 500,  negDashed = TRUE, curve = 0.2, posCol = "#808080",curveAll = TRUE, overlay = FALSE,  palette = "pastel", groups = phylum.group, vsize = 4, GLratio = 6, legend.cex = 0.35, overlay = TRUE)

##########################################
#RESISTANCE GENE NETWORK ANALYSIS (-rplB)#
##########################################

#perform network analysis without rplB!
#remove column based on pattern (rplB)
otu_table_norm.slim.t.genes <- otu_table_norm.slim.t[, -grep("rplB", colnames(otu_table_norm.slim.t))]

#find correlations between contigs!
corr.genes.slim <- corr.test(otu_table_norm.slim.t.genes, method = "spearman", adjust = "fdr", alpha = 0.03)

## prepare network graphics
#Create list of shapes and colors for network
r.slim <- data.frame(corr.genes.slim$r)
r.slim$gene <- rownames(r.slim)
r.slim$gene <- gsub("arsC_", "arsC", r.slim$gene)
r.slim <- r.slim %>%
  separate(gene, c("Gene", "Number"), by = "_", remove = FALSE) %>%
  left_join(gene, by = "Gene") %>%
  select(Gene, Number, Group, gene, gene.color) %>%
  rename(OTU = gene) %>%
  mutate(Shape = "circle")

#label r appropriately
r.slim$Group[r.slim$Gene == "rplB"] <- "Organism"
r.slim$Shape[r.slim$Group == "Organism"] <- "square"
r.slim$Group[is.na(r.slim$Group)] <- "Metadata"
r.slim$Shape[r.slim$Group == "Metadata"] <- "diamond"
r.slim$Shape[r.slim$Group == "AntibioticResistance"] <- "triangle"

#make vector of shapes for network
shapes.slim <- as.vector(r.slim$Shape)

#make rplB and metadata gene colors
r.slim$gene.color[r.slim$Gene == "rplB"] <- "#fffac8"
r.slim$gene.color[r.slim$Group == "Metadata"] <- "#D3D3D3"

#add leading zeros to r otus
r.final.slim <- r.slim %>%
  separate(OTU, into = c("beginning", "number"), sep = "_")

r.final.slim$number <- sprintf("%04s", r.final.slim$number)
r.final.slim$OTU <- paste(r.final.slim$beginning, r.final.slim$number, sep = "_")

#join aesthetic (shape, taxonomy) information
aesthetics.slim <- r.final.slim %>%
  left_join(blast.final, by = c("OTU", "Gene")) %>%
  unique()

#add phylum colors to aesthetics
aesthetics.slim <- aesthetics.slim %>%
  left_join(phylum.colors, by = "phylum")

#make taxanomic groups for network
phylum.group.slim <- split(as.numeric(rownames(aesthetics.slim)), list(aesthetics.slim$phylum)) 
class.group.slim <- split(as.numeric(rownames(aesthetics.slim)), list(aesthetics.slim$class), drop=TRUE) 
genus.group.slim <- split(as.numeric(rownames(aesthetics.slim)), list(aesthetics.slim$genus), drop=TRUE) 

#make network of correlations
clust.network.slim.taxonomy <- qgraph(corr.genes.slim$r, minimum = "sig", sampleSize=13, 
                                      details = TRUE, layout = "spring",
                                      graph = "cor",label.cex = 0.5,
                                      alpha = 0.0001, graph = "fdr", labels = aesthetics.slim$Gene,  label.scale.equal = TRUE, label.scale = FALSE,shape = aesthetics.slim$Shape, node.resolution = 500,  negDashed = TRUE, curve = 0.2, posCol = "#808080",curveAll = TRUE, overlay = FALSE,  palette = "ggplot2", groups = phylum.group.slim, vsize = 4, GLratio = 6, legend.cex = 0.35)

clust.network.slim.gene <- qgraph(corr.genes.slim$r, minimum = "sig", sampleSize=13, 
                                  details = TRUE, layout = "spring",
                                  graph = "cor",label.cex = 0.5,
                                  alpha = 0.01, graph = "fdr", labels = aesthetics.slim$number,  label.scale.equal = TRUE, label.scale = FALSE,shape = aesthetics.slim$Shape, node.resolution = 500,  negDashed = TRUE, curve = 0.2, posCol = "#808080",curveAll = TRUE, overlay = FALSE,  color = aesthetics.slim$gene.color, vsize = 4, GLratio = 6, legend.cex = 0.35)
