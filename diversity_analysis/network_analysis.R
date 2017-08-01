library(vegan)
library(psych)
library(tidyverse)
library(qgraph)
library(phyloseq)
library(reshape2)

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

#read in rplB data
#rplB <- read_delim(paste(wd, "/output/rplB.summary.scg.txt", sep = ""), delim  = " ")

#get mean rplB
#mean.rplB <- mean(rplB$rplB)

#get rplB otu data
rplB <- otu_table[,grepl("rplB", names(otu_table))]

#add site names to rplB
rownames(rplB) <-  otu_table$Site

#summarise rplB by adding all OTU counts 
# in each row (total rplB/site)
rplB_summary <- data.frame(rowSums(rplB))

#make site a column in rplB
rplB_summary$Site <- rownames(rplB_summary) 


#add rplB data to otu_table
otu_table.rplB <- rplB_summary %>%
  left_join(otu_table, by = "Site") %>%
  rename(rplB = rowSums.rplB.)

#normalize to rplB
otu_table_norm <- otu_table.rplB
for(i in 3:6570){otu_table_norm[,i]=otu_table.rplB[,i]/otu_table.rplB[,1]}

#add in metadata
otu_table_norm_annotated <- otu_table_norm %>%
  left_join(meta, by = "Site") %>%
  mutate(DateSince_Fire = 2014-DateFire_Elick2011) %>%
  select(rplB:As_ppm, SoilTemperature_to10cm, DateSince_Fire,
         OrganicMatter_500:Fe_ppm)


#change to df and add row names back
otu_table_norm_annotated <- as.data.frame(otu_table_norm_annotated)
rownames(otu_table_norm_annotated) <- otu_table_norm_annotated[,2]

#remove first two columns
otu_table_norm_annotated=otu_table_norm_annotated[,-c(1,2)]

#make data matrix
otu_table_norm_annotated=data.matrix(otu_table_norm_annotated)

#transpose data
otu_table_norm_annotated.t <- t(otu_table_norm_annotated)

#replace NAs with zeros
otu_table_norm_annotated.t[1:6569,][is.na(otu_table_norm_annotated.t[1:6569,])] <- 0


#make presence absence matrix
otu_table_normPA <- (otu_table_norm_annotated.t>0)*1

#replace NA with 0 (date since fire)
otu_table_normPA[is.na(otu_table_normPA)] <- 0
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

##make network with few genes
#subset data to only contain all of the AsRG ARG clusters
#list otus that are gene matches
matches <- c("tolC_05", "ClassB_003", "dfra12_076", "dfra12_038",
             "acr3_002", "acr3_053", "arsM_109", "arsM_296", "rplB_1027", 
             "rplB_1015", "rplB_0549", "rplB_0564", "rplB_0692", 
             "rplB_0149", "rplB_0407", "rplB_0355", "rplB_1131", 
             "rplB_0267", "rplB_0169", "rplB_0886", "SoilTemperature_to10cm",
             "Ca_ppm", "Mg_ppm", "DateSince_Fire")

#remove OTUs that are gene matches
otu_table_norm.slim.t_matches <- otu_table_norm.slim.t[,colnames(otu_table_norm.slim.t) %in% matches]

#find correlations between OTUs
corr <- corr.test(otu_table_norm.slim.t_matches, 
                  method = "spearman", adjust = "fdr", alpha = 0.01)

corr.r <- as.matrix(print(corr$r, long = TRUE))
corr.p <- as.matrix(print(corr$p, long = TRUE))

corr.r[which(corr.p > 0.01)] <- 0

#save correlation data
write.table(corr$r, paste(wd, "/output/corr_table.rplB.txt", sep = ""), quote = FALSE)

#make network of correlations
qgraph(corr$r, minimum = "sig", sampleSize=13, 
       layout = "circle", details = TRUE,
       graph = "cor", label.cex = 1,
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
gene_abundance <- gene_abundance[-which(gene_abundance$OTU == "rplB"),]

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

#order based on group
gene_abundance_summary$Gene <- factor(gene_abundance_summary$Gene, 
                                levels = gene_abundance_summary$Gene[order(gene_abundance_summary$Description)])

(boxplot <- ggplot(subset(gene_abundance_summary, subset = Group == "ArsenicResistance"), aes(x = Classification, 
                                               y = Total*100)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm)) +
    facet_wrap(~Gene, scales = "free_y", ncol = 3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                         guide_legend(title="Temperature (°C)")) +
    ylab("Gene per rplB (%)") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

ggsave(boxplot, filename = "/Users/dunivint/Documents/ShadeLab/Presentations/ARG_symposium/gene_boxplots_asrg.eps", height = 6, width = 7)

#subset data to only contain all of the AsRG ARG clusters
#list otus that are gene matches
matches <- c("tolC_05", "ClassB_003", "dfra12_076", "dfra12_038",
             "acr3_002", "acr3_053", "arsM_109", "arsM_296", "rplB_1027", 
             "rplB_1015", "rplB_0549", "rplB_0564", "rplB_0692", 
             "rplB_0149", "rplB_0407", "rplB_0355")

#remove OTUs that are gene matches
gene_abundance_mixes <- gene_abundance[gene_abundance$OTU %in% matches,]

#order based on temperature
gene_abundance_mixes$Site <- factor(gene_abundance_mixes$Site, 
                                      levels = gene_abundance_mixes$Site[order(gene_abundance_mixes$RelativeAbundance)])

ggplot(gene_abundance_mixes, aes(x = SoilTemperature_to10cm, 
                                 y = RelativeAbundance*100, fill = OTU)) +
  geom_density(stat = "identity", alpha = 0.5, position = "stack") +
  geom_point(aes(shape = Classification)) +
  facet_wrap(~OTU) +
  theme_bw()

#subset data to only contain all of the AsRG ARG clusters
#list otus that are gene matches
asrg <- c("acr3_160", "acr3_094", "acr3_047", "acr3_101", "acr3_018", 
             "acr3_024", "acr3_048", "arsM_023", "arsM_465", "arsM_101",
             "arsM_059", "arsM_072", "arsM_073", "arsM_262", "arsM_120",
             "arsCglut_094", "arsCglut_185", "arsCglut_066", "arsCglut_055",
             "arsCglut_139", "arsA_014", "arsA_010")
asrg <- c("acr3_022", "arsM_018", "arsCthio_06")
#remove OTUs that are gene matches
gene_abundance_asrg <- gene_abundance[gene_abundance$OTU %in% asrg,]

#order based on temperature
gene_abundance_asrg$Site <- factor(gene_abundance_asrg$Site, 
                                    levels = gene_abundance_asrg$Site[order(gene_abundance_asrg$RelativeAbundance)])

ggplot(gene_abundance_asrg, aes(x = Site, 
                                 y = RelativeAbundance*100, fill = OTU)) +
  geom_density(stat = "identity", alpha = 0.5) +
  geom_point(aes(shape = Classification)) +
  facet_wrap(~OTU) +
  theme_bw()


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

##############################
#TAXANOMIC ABUNDANCE ANALYSIS#
##############################

#need to adjust OTU column in gene_abundance
#to join with taxanomic data
gene_abundance_otu <- gene_abundance %>%
  separate(col = OTU, into = c("leftover", "OTU"), sep = "_") %>%
  select(-leftover)

#make OTU a number
gene_abundance_otu$OTU <- as.numeric(gene_abundance_otu$OTU)

#add leading zero to 4 digits
gene_abundance_otu$OTU <- sprintf("%04d", gene_abundance_otu$OTU)

#add OTU to otu label
gene_abundance_otu$OTU <- paste("OTU_", gene_abundance_otu$OTU, sep="")

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
blast.names <- list.files(pattern="blast_*")
identifiers <- do.call(rbind, lapply(blast.names, function(X) {
  data.frame(id = basename(X), read_delim(X, delim = "\t", col_types = 
                                            list(col_character(), 
                                                 col_character(), 
                                                 col_character(), 
                                                 col_character()),
                                          col_names = c("OTU", "Description",
                                                        "Evalue")))}))

#move back up a directory to proceed with analysis
setwd("../")

#Remove excess information in id column
identifiers$id <- gsub("blast_", "", identifiers$id)
identifiers$id <- gsub(".txt", "", identifiers$id)

#Tidy identifier data
identifiers_tidy <- identifiers %>%
  separate(col = Description, into = c("Accno", "Description"), 
           sep = " coded_by=") %>%
  separate(col = Description, into = c("Coded_by", "Description"),
           sep = ",organism=") %>%
  separate(col = Description, into = c("Taxon", "Definition"), 
           sep = ",definition=") %>%
  rename(Gene = id) %>%
  select(Gene, OTU, Taxon)

#list otus that are gene matches
mixed.aioA <- c("OTU_0013", "OTU_0031", "OTU_0034", "OTU_0036", "OTU_0037", "OTU_0042", "OTU_0043", "OTU_0049", "OTU_0051")
mixed.arrA <- c("OTU_0001", "OTU_0002", "OTU_0003", "OTU_0004", "OTU_0005", "OTU_0006", "OTU_0008")
mixed.arxA <- c("OTU_0011", "OTU_04")

#remove OTUs that are gene matches
identifiers_tidy <- identifiers_tidy[-which(identifiers_tidy$Gene == "aioA" &
                                              identifiers_tidy$OTU %in% mixed.aioA),]
identifiers_tidy <- identifiers_tidy[-which(identifiers_tidy$Gene == "arrA" &
                                              identifiers_tidy$OTU %in% mixed.arrA),]
identifiers_tidy <- identifiers_tidy[-which(identifiers_tidy$Gene == "arxA" &
                                              identifiers_tidy$OTU %in% mixed.arxA),]

library(taxize)
#add taxanomic information 
#blast.ncbi <- tax_name(query = identifiers_tidy$Taxon, 
#                      get = c("genus", "order", "family", "class", "phylum"), #db = "ncbi")


#label query "Organism" for joining purposes
#blast.ncbi$Taxon <- blast.ncbi$query

#save this table since the above step takes a long time
#write.table(blast.ncbi, file = paste(wd, "/output/blast.ncbi.taxonomy.txt",
#sep = ""), row.names = FALSE)

#read in ncbi information 
blast.ncbi <- read_delim(paste(wd, "/output/blast.ncbi.taxonomy.txt",
                               sep = ""), delim = " ")
  
#join ncbi information with annotated data
#output should have same number of rows 
identifiers_ncbi <- identifiers_tidy %>%
  left_join(blast.ncbi, by = "Taxon") %>%
  unique() %>%
  left_join(gene_abundance_otu, by = c("OTU", "Gene")) %>%
  unique()

#replace NA in phylum with unknown
identifiers_ncbi$phylum[is.na(identifiers_ncbi$phylum)] = "Unknown"

#call NA class by phyla
identifiers_ncbi$class[is.na(identifiers_ncbi$class)] <- as.character(identifiers_ncbi$phylum[is.na(identifiers_ncbi$class)])

#call NA genus by class (may be phyla in cases where class was NA)
identifiers_ncbi$genus[is.na(identifiers_ncbi$genus)] <- as.character(identifiers_ncbi$class[is.na(identifiers_ncbi$genus)])

#order based on temperature
identifiers_ncbi$Site <- factor(identifiers_ncbi$Site, 
                                   levels = identifiers_ncbi$Site[order(identifiers_ncbi$SoilTemperature_to10cm)])


##############################
#EXAMINE PHYLUM LEVEL CHANGES#
##############################

#look at phylum level differneces
data.phylum <- identifiers_ncbi %>%
  rename(Temp = SoilTemperature_to10cm) %>%
  group_by(Group, Description, Gene, phylum, Classification, Site, Temp) %>%
  summarise(Phylum.count = sum(RelativeAbundance))

#prep colors for phylum diversity
color <- c("#FF7F00", "#7570B3", "#CAB2D6", "#FBB4AE", "#F0027F", "#BEBADA", "#E78AC3", "#A6D854", "#B3B3B3", "#386CB0", "#BC80BD", "#FFFFCC", "#BF5B17", "#984EA3", "#CCCCCC", "#FFFF99", "#B15928", "#F781BF", "#FDC086", "#A6CEE3", "#FDB462", "#FED9A6", "#E6AB02", "#E31A1C", "#B2DF8A", "#377EB8", "#FCCDE5", "#80B1D3", "#FFD92F", "#33A02C", "#66C2A5", "#666666", "black", "brown", "grey", "red", "blue", "green", "purple", "brightorange", "pink", 
           "yellow")

#order genes by group
data.phylum$Gene <- factor(data.phylum$Gene, 
                           levels = data.phylum$Gene[order(data.phylum$Group)])

#plot 
(gene.bar.census <- ggplot(subset(data.phylum, Group == "ArsenicResistance"),
                           aes(x = Site,  y = Phylum.count*100, fill = phylum)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = color) +
    theme_classic(base_size = 8) +
    ylab("Gene per rplB (%)") +
    facet_wrap(~ Gene, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=0.95,vjust=0.2)))





