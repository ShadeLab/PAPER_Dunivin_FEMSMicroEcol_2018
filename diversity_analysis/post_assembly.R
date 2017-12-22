#######################################
#This script analyzes data assembled from Centralia
#metagenomes using Xander. It was written by
#Taylor K Dunivin

#For more information and input files see the 
#corresponding GitHub repo:
#https://github.com/ShadeLab/PAPER_Dunivin_Antibiotics_2017
#######################################

library(vegan)
library(psych)
library(tidyverse)
library(phyloseq)
library(reshape2)
library(broom)
library(data.table)

#print working directory for future references
wd <- print(getwd())

#setwd to diversity analysis
setwd(paste(wd, "/diversity_analysis", sep = ""))

#now change wd obj
wd <- print(getwd())

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_FULL_map.txt", sep=""), sep=" ", header=TRUE))

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
filenames <- list.files(pattern="*_rformat_dist_0.01.txt")

#move back up directories
setwd(wd)

#make dataframes of all OTU tables
for(i in filenames){
  filepath <- file.path(paste(wd, "/data", sep = ""),paste(i,sep=""))
  assign(gsub("_rformat_dist_0.01.txt", "", i), read.delim(filepath,sep = "\t"))
}

#change OTU to gene name
colnames(`AAC6-Ia`) <- naming(`AAC6-Ia`)
colnames(adeB) <- naming(adeB)
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
otu_table <- rplB %>%
  left_join(`AAC6-Ia`, by = "X") %>%
  left_join(CEP, by = "X") %>%
  left_join(ClassA, by = "X") %>%
  left_join(ClassB, by = "X") %>%
  left_join(ClassC, by = "X") %>%
  left_join(dfra12, by = "X") %>%
  left_join(intI, by = "X") %>%
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
  rename(Site = X) 

#replace all NAs (from join) with zeros
otu_table[is.na(otu_table)] <- 0

#write file to save OTU table
write.table(otu_table, paste(wd, "/output/otu_table.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

################################################
#EXTRACT AND TIDY RPLB FOR NORMALIZATION (rplB)#
################################################

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
write.table(rplB_summary_save, file = paste(wd, "/output/rplB.summary.scg.txt", sep = ""), row.names = FALSE)

#add rplB data to otu_table
otu_table.rplB <- rplB_summary %>%
  left_join(otu_table, by = "Site") %>%
  rename(rplB = rowSums.rplB.)

#normalize to rplB
otu_table_norm <- otu_table.rplB
for(i in 3:ncol(otu_table_norm)){otu_table_norm[,i]=otu_table.rplB[,i]/otu_table.rplB[,1]}

#remove Cen13, which is a chemical and community (16S outlier)
otu_table_norm <- otu_table_norm[!otu_table_norm$Site == "cen13",]

#add in metadata
otu_table_norm$Site <- gsub("cen", "Cen", otu_table_norm$Site)
otu_table_norm_annotated <- otu_table_norm %>%
  left_join(meta, by = "Site") %>%
  mutate(DateSince_Fire = 2014-DateFire_Elick2011) %>%
  select(rplB:As_ppm, SoilTemperature_to10cm, DateSince_Fire,
         OrganicMatter_500:SoilMoisture_perc)

#change to df and add row names back
#remove first two columns (redundant)
otu_table_norm_annotated <- as.data.frame(otu_table_norm_annotated)
rownames(otu_table_norm_annotated) <- otu_table_norm_annotated[,2]
otu_table_norm_annotated=otu_table_norm_annotated[,-c(1,2)]

#make data matrix and transpose 
otu_table_norm_annotated <- data.matrix(otu_table_norm_annotated)
otu_table_norm_annotated.t <- t(otu_table_norm_annotated)

#replace NAs with zeros (except date since fire)
n <- nrow(otu_table_norm_annotated.t)-12
otu_table_norm_annotated.t[1:n,][is.na(otu_table_norm_annotated.t[1:n,])] <- 0

#make presence absence matrix
otu_table_normPA <- (otu_table_norm_annotated.t>0)*1

#replace NA with 0 (date since fire)
otu_table_normPA[is.na(otu_table_normPA)] <- 0

########################
#ARG TAXANOMIC ANALYSIS#
########################

#list OTUs present in 2 or more samples
abund_3 <- otu_table_normPA[which(rowSums(otu_table_normPA) >2),]

#remove OTUs with presence in 2 or less samples
otu_table_norm.slim.t_3 <- t(otu_table_norm_annotated.t[which(rownames(otu_table_norm_annotated.t) %in% rownames(abund_3)),])
otu_table_norm.slim.t <- t(otu_table_norm_annotated.t)

#tidy aesthetics and prepare for blast join
aesthetics_3 <- otu_table_norm.slim.t_3 %>%
  melt(value.name = "normAbund") %>%
  rename(Site = Var1, OTU = Var2) %>%
  left_join(meta, by = "Site") %>%
  separate(OTU, into = c("beginning", "number"), sep = "_") 
aesthetics_3$number <- sprintf("%04s", aesthetics_3$number)
aesthetics_3$OTU <- paste(aesthetics_3$beginning, aesthetics_3$number, sep = "_")

aesthetics <- otu_table_norm.slim.t %>%
  melt(value.name = "normAbund") %>%
  rename(Site = Var1, OTU = Var2) %>%
  left_join(meta, by = "Site") %>%
  separate(OTU, into = c("beginning", "number"), sep = "_") 
aesthetics$number <- sprintf("%04s", aesthetics$number)
aesthetics$OTU <- paste(aesthetics$beginning, aesthetics$number, sep = "_")

#add taxanomic information to OTUs
#temporarily change working directories
#setwd(paste(wd, "/data/BLAST_results", sep = ""))

#list filenames of interest
#blast.filenames <- list.files(pattern="*0.01.txt")

#read in blast data
#blast.data <- do.call(rbind, lapply(blast.filenames, function(X) {  data.frame(Gene = gsub("results_|_0.01.txt|_", "",basename(X)), read_delim(X, delim = "\t", col_names = FALSE))}))

#move back up directories
#setwd("../../")

#separate out rplB from other genes
#blast.rplB <- blast.data[blast.data$Gene == "rplB",]
#blast.all <- blast.data[!blast.data$Gene == "rplB",]

#tidy blast data
#blast.all$X1 <- gsub("OTU_", "", blast.all$X1)
#blast.all.tidy <- blast.all %>%  mutate(OTU = paste(Gene, X1, sep = "_")) %>%  separate(X3, into = c("Accno", "Extra"), sep = " coded_by=") %>%  separate(Extra, into = c("Complement", "Organism"), sep = ",organism=") %>%  separate(Organism, into = c("Organism", "Description"), sep = ",definition=") 

#remove all rows with uncultured orgs
#blast.all.tidy <- blast.all.tidy[!grepl("uncultured|metagenome", blast.data.tidy$Organism),]
#blast.all.tidy <- blast.all.tidy[!grepl("arsA|arxA|arsM|aioA|arsC|arrA|arsB|acr3", blast.all.tidy$Gene),]
#blast.unique <- unique(blast.all.tidy$Organism)

#get taxanomic information for blasted organisms
#library(taxize)
#taxize_results <- tax_name(query = blast.unique, get = c("phylum", "class", "genus"), db = "ncbi")
#taxize_results$phylum[grep("Parcubacteria", taxize_results$query)] <- "Parcubacteria"
#taxize_results$phylum[grep("Microgenomates", taxize_results$query)] <- "Microgenomates"
#taxize_results$phylum[grep("Cupriavidus", taxize_results$query)] <- "Proteobacteria"
#taxize_results$class[grep("Cupriavidus", taxize_results$query)] <- "Betaproteobacteria"
#taxize_results$genus[grep("Cupriavidus", taxize_results$query)] <- "Cupriavidus"
#taxize_results$phylum[grep("aminovorans", taxize_results$query)] <- "Firmicutes"
#taxize_results$class[grep("aminovorans", taxize_results$query)] <- "Bacilli"
#taxize_results$genus[grep("aminovorans", taxize_results$query)] <- "aminovorans"
#taxize_results$phylum[grep("X1", taxize_results$query)] <- "Firmicutes"
#taxize_results$phylum[grep("Polyangium", taxize_results$query)] <- "Proteobacteria"
#taxize_results$class[grep("Polyangium", taxize_results$query)] <- "Deltaproteobacteria"
#taxize_results$genus[grep("Polyangium", taxize_results$query)] <- "Polyangium"
#taxize_results$phylum[grep("PML1", taxize_results$query)] <- "Proteobacteria"
#taxize_results$class[grep("PML1", taxize_results$query)] <- "Betaproteobacteria"
#taxize_results$genus[grep("PML1", taxize_results$query)] <- "Burkholderia"
#taxize_results$class[is.na(taxize_results$class)] <- taxize_results$phylum[is.na(taxize_results$class)]
#taxize_results$phylum <- ifelse((taxize_results$phylum == "Proteobacteria"), taxize_results$class, taxize_results$phylum)

#blast.all.tidy.tax <- blast.all.tidy %>%  rename(query = Organism) %>%  left_join(taxize_results, by = "query") 
#blast.all.tidy.tax <- blast.all.tidy.tax[!duplicated(blast.all.tidy.tax$OTU),] %>%  select(Gene, OTU, phylum, class, genus)

#set up rplB 
#blast.rplB.tidy <- blast.rplB %>%  separate(X3, into = c("organism", "phylum", "class", "order", "family", "genus"), sep = ";") %>%  select(Gene, X1, phylum, class, genus) %>% rename(OTU = X1)
#blast.rplB.tidy$OTU <- gsub("OTU_", "rplB_", blast.rplB.tidy$OTU)

#combine gene and rplB data
#blast.final <- rbind(blast.all.tidy.tax, blast.rplB.tidy)

#save results to file
#write.table(blast.final, paste(wd, "/output/taxize_0.01_results_FINAL.txt", sep = ""), col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

#read in taxize / blast information
blast.final <- read_delim(paste(wd, "/output/taxize_0.01_results_FINAL.txt", sep = ""), col_names = TRUE, delim = "\t")

#join blast data with aesthetics
#for both "abundant" OTUs and all OTUs
aesthetics.blast_3 <- aesthetics_3 %>%
  left_join(blast.final, by = "OTU")
aesthetics.blast_3 <- aesthetics.blast_3[!duplicated(aesthetics.blast_3[c("OTU","Site")]),]

#read in gene classification data
gene <- read_delim(paste(wd, "/data/gene_classification.txt",  sep=""), 
                   delim = "\t", col_names = TRUE)

aesthetics.blast <- aesthetics %>%
  left_join(blast.final, by = "OTU") %>%
  left_join(gene, by = "Gene")
aesthetics.blast <- aesthetics.blast[!duplicated(aesthetics.blast[c("OTU","Site")]),]

#test for changes in abundant OTUs along chronosequence
#against temperature
spearman_OTUtemp <- subset(aesthetics.blast_3, Gene !="rplB") %>% group_by(Gene, OTU) %>% do(tidy(cor.test(.$normAbund, log(.$SoilTemperature_to10cm), method = "spearman")))
write.table(spearman_OTUtemp, file = paste(wd, "/output/spearman_OTUtemp.csv", sep = ""), quote = FALSE, row.names = FALSE, sep = ",")

#summarise and all OTUs based on class
aesthetics.blast.summary <- aesthetics.blast %>%
  group_by(Site, SoilTemperature_to10cm, Gene, phylum) %>%
  summarise(PhyTot = sum(normAbund)) %>%
  ungroup() %>%
  group_by(Site, Gene) %>%
  mutate(SiteTot = sum(PhyTot), RelAbund = PhyTot/SiteTot)

#test for changes in ARG phylum (taxonomy) along chronosequence
#against temperature
spearman_Phytemp <- subset(aesthetics.blast.summary) %>% group_by(Gene, phylum) %>% do(tidy(cor.test(.$PhyTot, log(.$SoilTemperature_to10cm), method = "spearman")))
write.table(spearman_Phytemp, file = paste(wd, "/output/spearman_Phytemp.csv", sep = ""), quote = FALSE, row.names = FALSE, sep = ",")

#test for changes in ARG NORMALIZED phylum (taxonomy) along #chronosequence against temperature
spearman_nPhytemp <- subset(aesthetics.blast.summary) %>% group_by(Gene, phylum) %>% do(tidy(cor.test(.$RelAbund, log(.$SoilTemperature_to10cm), method = "spearman")))
write.table(spearman_nPhytemp, file = paste(wd, "/output/spearman_nPhytemp.csv", sep = ""), quote = FALSE, row.names = FALSE, sep = ",")

#examine changes in class level abundance 
list <- c("ClassA", "ClassB","dfra12", "intI", "rplB")

#change proteobacteria to class
aesthetics.blast$phylum <- ifelse((aesthetics.blast$phylum == "Proteobacteria"), aesthetics.blast$class, aesthetics.blast$phylum)
aesthetics.blast$phylum <- ifelse((aesthetics.blast$phylum == "environmentalsamples"), "metagenomes", aesthetics.blast$phylum)
aesthetics.blast$phylum <- ifelse((aesthetics.blast$phylum == "artificialsequences"), "metagenomes", aesthetics.blast$phylum)

#order by temperature
aesthetics.blast$Site <- factor(aesthetics.blast$Site, 
                                aesthetics.blast$Site[order(aesthetics.blast$SoilTemperature_to10cm)])

colors37 = c("#466791","#953ada","#4fbe6c","#ce49d3","#a7b43d","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975","grey90")

#summarise dataset
summary.aesthetics.blast <- aesthetics.blast %>%
  group_by(Site, SoilTemperature_to10cm, Gene, phylum) %>%
  summarise(sumNormAbund = sum(normAbund))
(class.prop <- ggplot(subset(summary.aesthetics.blast, Gene %in% list), aes(x = Site, y = sumNormAbund)) +
    geom_bar(stat = "identity", position = "fill", aes(fill = phylum)) +
    facet_wrap(~Gene, ncol = 2) +
    scale_fill_manual(values = colors37) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=0.95,vjust=0.9)))
ggsave(class.prop, filename = paste(wd, "/figures/class_proportions_temp.eps", sep = ""), width = 12, height= 6.5,units = "in")

###############################
#MAKE HEATMAP TO TEST PATTERNS#
###############################

#list OTUs present in 2 or more samples
abund_2 <- otu_table_normPA[which(rowSums(otu_table_normPA) >1),]

#remove OTUs with presence in 2 or less samples
otu_table_norm.slim.t_2 <- data.frame(t(otu_table_norm_annotated.t[which(rownames(otu_table_norm_annotated.t) %in% rownames(abund_2)),]))
otu_table_norm.slim.t_2$Site <- rownames(otu_table_norm.slim.t_2)

#order dataframe
otu_table_norm.slim.t_2 <- arrange(otu_table_norm.slim.t_2,SoilTemperature_to10cm)
rownames(otu_table_norm.slim.t_2) <- otu_table_norm.slim.t_2$Site
otu_table_norm.slim.t_2 <- as.matrix(otu_table_norm.slim.t_2[,c(1:349)])
otu_table_norm.slim.t_2 <- otu_table_norm.slim.t_2[,-grep("rplB", colnames(otu_table_norm.slim.t_2))]

#get heatmap colors
hc=colorRampPalette(c("white", "#91bfdb", "midnightblue"), interpolate="linear", bias = 3)

#prep data for gene annotation on heatmap
colors.otu.2 <- data.frame(t(otu_table_norm.slim.t_2))
colors.otu.2_annotated <- colors.otu.2 %>%
  rownames_to_column(var = "OTU") %>%
  separate(col = OTU, into = c("Gene", "OTU"), sep = "_") %>%
  mutate(Max = apply(colors.otu.2, 1, max)) %>%
  select(Gene)
rownames(colors.otu.2_annotated) <- rownames(colors.otu.2)


#set up environment to run heatmap
library(pheatmap)
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

#set gene colors for plotting
ann_colors = list(
  Gene = c(adeB = "#8DD3C7", CEP = "#FFFFB3", ClassA = "#BEBADA", ClassB = "#FB8072", ClassC = "#80B1D3", dfra12 = "#FDB462", intI = "#B3DE69", sul2 = "#FCCDE5", tolC = "#D9D9D9", vanA = "#BC80BD", vanH = "#CCEBC5", vanX = "chocolate4", vanZ = "grey40"))

#plot heatmap
(heatmap <- pheatmap(t(otu_table_norm.slim.t_2), cluster_rows = TRUE, cluster_cols = FALSE, clustering_method = "complete", dendrogram = "row", scale = "none", trace = "none", legend = TRUE, color = hc(5000), cellheight = 9, cellwidth = 24, treeheight_row = 250, fontsize = 16, border_color = NA, show_rownames = FALSE, annotation_row = colors.otu.2_annotated, annotation_colors = ann_colors, clustering_callback = callback, width = 12, height = 16))

#heatmaps/dendrograms show several groups (listed below)

######################
#CORRELATION ANALYSES#
######################
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

#order based on group
gene_abundance_summary$Gene <- factor(gene_abundance_summary$Gene, 
                                      levels = gene_abundance_summary$Gene[order(gene_abundance_summary$Description)])

#decast for abundance check but first make Site a character
gene_abundance_summary$Site <- as.character(gene_abundance_summary$Site)
cast.gene <- data.frame(acast(gene_abundance_summary, Site ~ Gene, value.var = "Total"))

#call na's zeros
cast.gene[is.na(cast.gene)] =0

#remove unnecessary sites in meta data
meta.slim <- meta[meta$Site %in% gene_abundance_summary$Site,]

#correlate data with temp
corr <- print(corr.test(x = cast.gene, y = meta.slim[,c(2,16, 20:30,32)], method = "spearman", adjust = "fdr"), short = FALSE)
cor <- corr.test(x = cast.gene, y = meta.slim[,c(2,16, 20:30,32)], method = "spearman", adjust = "fdr")
#antibiotic resistance genes are not correlated with env variables!

#save as table for supplemental material
write.table(corr, file = paste(wd, "/output/geochem.correlations.csv", sep = ""), sep = ",", row.names = TRUE, quote = FALSE)

#correlate genes with log(Temp)
spearman_logtemp <- gene_abundance_summary %>% group_by(Gene) %>% do(tidy(cor.test(.$Total, log(.$SoilTemperature_to10cm), method = "spearman")))

#save as table for supplemental material
write.table(spearman_logtemp, file = paste(wd, "/output/logTemp.correlations.csv", sep = ""), sep = ",", row.names = FALSE, quote = FALSE)

#check correlations between genes
cast.gene <- cast.gene %>% select(-c(rplB,tetA, tetW, tetX))
cast.gene <- cast.gene[,-1]
corr.genes <- print(corr.test(cast.gene, method = "spearman", adjust = "fdr"), 
                    short = FALSE)
corr.genes.matrix <- corr.test(cast.gene, method = "spearman", adjust = "fdr")

#save correlations to table
write.table(corr.genes, paste(wd, "/output/gene.correlations.csv", sep = ""), row.names = TRUE, sep = ",", quote = FALSE)

cor.plot(corr.genes.matrix$r,numbers = TRUE, xlas = 2, upper = FALSE, diag = FALSE, stars = TRUE, pval = corr.genes.matrix$p)

############################
#MAKE GENE ABUNDANCE GRAPHS#
############################

sig <- c("ClassA", "ClassB", "dfra12", "tolC")
#plot antibiotic resistance genes
(sigTempGene <- ggplot(subset(gene_abundance_summary, subset = Gene %in% sig), aes(x = SoilTemperature_to10cm, 
                                                                                   y = Total)) +
    geom_smooth(method = "lm", color = "grey50", se = FALSE) +
    geom_point(aes(shape = Classification), size = 3) +
    facet_wrap(~Gene, scales = "free_y", ncol = 2) +
    ylab("rplB-normalized abundance") +
    scale_x_log10("Soil Temperature", breaks = c(12,20,30,40,50)) +
    theme_bw(base_size = 16))


ggsave(sigTempGene, filename = paste(wd, "/figures/sig_temp_gene.eps", sep = ""), height = 5, width = 7, units = "in")

(nonTempGene <- ggplot(subset(gene_abundance_summary, subset = !Gene %in% sig), aes(x = SoilTemperature_to10cm, 
                                                                                    y = Total)) +
    geom_smooth(method = "lm", color = "grey50", se = FALSE) +
    geom_point(aes(shape = Classification), size = 3) +
    facet_wrap(~Gene, scales = "free_y", ncol = 3) +
    ylab("rplB-normalized abundance") +
    scale_x_log10("Soil Temperature (°C)", breaks = c(0,10,20,30,40,50)) +
    theme_bw(base_size = 12))

ggsave(nonTempGene, filename = paste(wd, "/figures/nonTempGene.eps", sep = ""), height = 8, width = 7, units = "in")

####################################
#EXAMINE rplB ACROSS CHRONOSEQUENCE#
####################################

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names <- list.files(pattern="*_45_taxonabund.txt")
data <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read_delim(X, delim = "\t"))}))

#move back up a directory to proceed with analysis
setwd(wd)

#split columns, tidy dataset, and combine proteobacteria
data <- data %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = junk, into = c("Gene", "junk"), sep = "_45_", remove = TRUE) %>%
  mutate(Gene = gsub("_", "", Gene),
         Site = gsub("cen", "Cen", Site)) 

#split columns and set abundance and fract. abundace
#to numbers instead of characters
rplB <- data %>%
  select(Site, Taxon:Fraction.Abundance) %>%
  group_by(Site) %>%
  mutate(Fraction.Abundance = as.numeric(Fraction.Abundance),
         Abundance = as.numeric(Abundance)) %>%
  mutate(Taxon = gsub(".*proteobacteria", "Proteobacteria", Taxon)) %>%
  subset(!Site == "Cen13")

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.rplB <- rplB %>%
  summarise(Total = sum(Fraction.Abundance), rplB = sum(Abundance))

#decast for abundance check and call na's zero
dcast <- acast(rplB, Taxon ~ Site, value.var = "Fraction.Abundance", fun.aggregate = sum)
dcast[is.na(dcast)] = 0

#order based on abundance
order.dcast <- dcast[order(rowSums(dcast),decreasing=TRUE),]

#melt data
melt <- melt(order.dcast,id.vars=row.names(order.dcast), variable.name= "Site", value.name="Fraction.Abundance" )

#adjust colnames of melt
colnames(melt) <- c("Taxon", "Site", "Fraction.Abundance")

#join metadata with regular data
history <- melt %>%
  left_join(meta, by="Site") %>%
  group_by(Taxon, Classification) %>%
  summarise(N = length(Fraction.Abundance), 
            Average = mean(Fraction.Abundance))

#plot
(phylum.plot=(ggplot(history, aes(x=Taxon, y=Average)) +
                geom_point(size=2) +
                facet_wrap(~Classification, ncol = 1) +
                labs(x="Phylum", y="Mean relative abundance")+
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, size = 10, 
                                                 hjust=0.95,vjust=0.2))))

#save plot
ggsave(phylum.plot, filename = paste(wd, "/figures/phylum.responses.eps", sep=""), 
       width = 5, height = 5, units = "in")

#test rplB correlations with temperature
#(norm abund)
spearman_rplBtemp <- melt %>% left_join(meta, by = "Site") %>% group_by(Taxon) %>% do(tidy(cor.test(.$Fraction.Abundance, log(.$SoilTemperature_to10cm), method = "spearman")))

#########################################################
#a-DIVERSITY ANALYSIS OF ARGs AND AsRGs clustered at 0.01#
#########################################################

#separate original (non-normalized) contig table into
#specific gene groups (rplB, ARG, AsRG)
#remove cen13!
otu_table <- otu_table[!otu_table$Site == "cen13",]
ecol.rplB <- otu_table[, grep("rplB", colnames(otu_table))]
ecol.ARG <- otu_table[,grep("tolC|dfra|Class|tet|van|CEP|AAC|ade|sul", colnames(otu_table))]

#add site names as row names to each contig table 
rownames(ecol.rplB) <- otu_table[,1]
rownames(ecol.ARG) <- otu_table[,1]

#make list of 13 colors (based on classification)
class.13 <- c("red", "yellow", "red", "yellow", "red", "red","darkgreen", "yellow", "yellow", "red", "yellow", "red")

#check sampling depth of each matrix
rarecurve(ecol.rplB, step=1, label = FALSE, col = class.13)
rarecurve(ecol.ARG, step=1, label = FALSE, col = class.13)

#convert contig tables into phyloseq objects
ecol.rplB.phyloseq <- otu_table(ecol.rplB, taxa_are_rows = FALSE)
ecol.ARG.phyloseq <- otu_table(ecol.ARG, taxa_are_rows = FALSE)

#rarefy rplB to even sampling depth 
ecol.rplB.rare <- rarefy_even_depth(ecol.rplB.phyloseq, rngseed = TRUE)
rarecurve(ecol.rplB.rare, step=1, label = FALSE, col = class.13)

#rarefy ARG to even sampling depth
ecol.ARG.rare <- rarefy_even_depth(ecol.ARG.phyloseq, rngseed = TRUE)
rarecurve(ecol.ARG.rare, step=1, label = FALSE, col = class.13)

#calculate evenness
plieou.rplB <- data.frame(group = "rplB", Site = rownames(ecol.rplB.rare), plieou = vegan::diversity(ecol.rplB.rare, index = "shannon")/log(specnumber(ecol.rplB.rare)))

plieou.ARG <- data.frame(group = "ARG", Site = rownames(ecol.ARG.rare), plieou = vegan::diversity(ecol.ARG.rare, index = "shannon")/log(specnumber(ecol.ARG.rare)))


#join all evenness information and add metadata
plieou.full <- rbind(plieou.ARG, plieou.rplB)
plieou.full <- plieou.full %>%
  mutate(Site = gsub("cen", "Cen", Site)) %>%
           left_join(meta, by = "Site")

#plot evenness
(plieou.plot <- ggplot(plieou.full, aes(x = SoilTemperature_to10cm, y = plieou)) +
    geom_smooth(method = "lm", se = FALSE, color = "grey50") +
    geom_point(aes(shape = Classification), size = 3) +
    ylab(label = "Evenness") +
    theme_bw(base_size = 20) +
    facet_wrap(~group))

#save evennes plots
ggsave(plieou.plot, filename = paste(wd, "/figures/evenness.eps", sep = ""), width = 7.5, units = "in", height = 3)

#correlate evenness with Temp
spearman_evenness <- plieou.full %>% group_by(group) %>% do(tidy(cor.test(.$plieou, .$SoilTemperature_to10cm, method = "spearman")))

#save as table for supplemental material
write.table(spearman_evenness, file = paste(wd, "/output/evenness.Temp.correlations.csv", sep = ""), sep = ",", row.names = FALSE, quote = FALSE)

#make metadata a phyloseq class object
meta$Site <- gsub("Cen", "cen", meta$Site)
rownames(meta) <- meta[,1]
meta.phylo <- meta[,-1]
meta.phylo <- sample_data(meta.phylo)

##make biom for phyloseq
ecol.rplB.rare <- merge_phyloseq(ecol.rplB.rare, meta.phylo)
ecol.ARG.rare <- merge_phyloseq(ecol.ARG.rare, meta.phylo)

#plot & save RICHNESS
(richness.ecol.rplB.rare <- plot_richness(ecol.rplB.rare, x = "SoilTemperature_to10cm", measures = "Observed") +
    geom_smooth(method = "lm", se = FALSE, color = "grey50") +
    geom_point(aes(shape = Classification), size = 3) +
    ylim(35,175) +
    theme_bw(base_size = 18))
ggsave(richness.ecol.rplB.rare, filename = paste(wd, "/figures/richness.rplB.eps", sep = ""), width = 5, height = 3, units = "in")

(richness.ecol.ARG.rare <- plot_richness(ecol.ARG.rare, x = "SoilTemperature_to10cm", measures = "Observed") +
    geom_smooth(method = "lm", se = FALSE, color = "grey50") +
    geom_point(aes(shape = Classification), size = 3) +
    ylim(35,175) +
    theme_bw(base_size = 18))
ggsave(richness.ecol.ARG.rare, filename = paste(wd, "/figures/richness.ARG.eps", sep = ""), width = 5, height = 3, units = "in")


#extract richness and perform statistical tests
richness.ecol.rplB.rare.data <- estimate_richness(ecol.rplB.rare, split = TRUE, measures = "Observed")
richness.ecol.rplB.rare.data <- richness.ecol.rplB.rare.data %>%
  mutate(Site = rownames(richness.ecol.rplB.rare.data), GeneGroup = "rplB")

richness.ecol.ARG.rare.data <- estimate_richness(ecol.ARG.rare, split = TRUE, measures = "Observed")
richness.ecol.ARG.rare.data <- richness.ecol.ARG.rare.data %>%
  mutate(Site = rownames(richness.ecol.ARG.rare.data), GeneGroup = "ARG")

#join together richness results and add metadata 
richness.data <- rbind(richness.ecol.rplB.rare.data, richness.ecol.ARG.rare.data)
richness.data <- richness.data %>%
  left_join(meta, by = "Site")

#correlate richness with Temp
spearman_richness <- richness.data %>% group_by(GeneGroup) %>% do(tidy(cor.test(.$Observed, .$SoilTemperature_to10cm, method = "spearman")))

#save as table for supplemental material
write.table(spearman_richness, file = paste(wd, "/output/richness.Temp.correlations.csv", sep = ""), sep = ",", row.names = FALSE, quote = FALSE)

#relativize rarefied datasets before beta-diversity
ecol.rplB.rareREL <-  transform_sample_counts(ecol.rplB.rare, function(x) x/sum(x))
ecol.ARG.rareREL <-  transform_sample_counts(ecol.ARG.rare, function(x) x/sum(x))

#########################################################
#B-DIVERSITY ANALYSIS OF ARGs AND AsRGs clustered at 0.01#
#########################################################

#plot Bray Curtis ordination for rplB
ord.rplB.bray <- ordinate(ecol.rplB.rareREL, method="PCoA", distance="bray")
(bc.ord.rplB=plot_ordination(ecol.rplB.rareREL, ord.rplB.bray, color="SoilTemperature_to10cm",
                             title="Bray Curtis - rplB", shape = "Classification") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    geom_point(size=5) +
    theme_light(base_size = 12))

#plot Bray Curtis ordination for ARGs
ord.ARG.bray <- ordinate(ecol.ARG.rareREL, method="PCoA", distance="bray")
(bc.ord.ARG=plot_ordination(ecol.ARG.rareREL, ord.ARG.bray, color="SoilTemperature_to10cm",
                            title="Bray Curtis - ARG", shape = "Classification") +
    geom_point(size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save each bray curtis
ggsave(bc.ord.rplB, filename = paste(wd, "/figures/rplB.braycurtis.eps", sep = ""), width = 5, height =3, units = "in")
ggsave(bc.ord.ARG, filename = paste(wd, "/figures/ARG.braycurtis.eps", sep = ""), width = 5, height = 3, units = "in")

#extract OTU table from phyloseq
ecol.rplB.rareREL.matrix = as(otu_table(ecol.rplB.rareREL), "matrix")
ecol.ARG.rareREL.matrix = as(otu_table(ecol.ARG.rareREL), "matrix")

#calculate distance matrix
ecol.rplB.rareREL.d <- vegdist(ecol.rplB.rareREL.matrix, diag = TRUE, upper = TRUE)
ecol.ARG.rareREL.d <- vegdist(ecol.ARG.rareREL.matrix, diag = TRUE, upper = TRUE)

#mantel tests
mantel(ecol.rplB.rareREL.d,ecol.ARG.rareREL.d, method = "spear")

#try mantel tests of recovered sites only
recovered <- c("cen01", "cen03", "cen04", "cen05", "cen07")

ecol.rplB.rareREL.Rec.d <- vegdist(ecol.rplB.rareREL.matrix[rownames(ecol.rplB.rareREL.matrix) %in% recovered,], diag = TRUE, upper = TRUE)

ecol.ARG.rareREL.Rec.d <- vegdist(ecol.ARG.rareREL.matrix[rownames(ecol.ARG.rareREL.matrix) %in% recovered,], diag = TRUE, upper = TRUE)

#mantel tests on recovered only
mantel(ecol.rplB.rareREL.Rec.d,ecol.ARG.rareREL.Rec.d, method = "spear")

#try mantel tests on fire affected sites only
fireaffected <- c("cen06", "cen10", "cen12", "cen14", "cen15", "cen16")

ecol.rplB.rareREL.FA.d <- vegdist(ecol.rplB.rareREL.matrix[rownames(ecol.rplB.rareREL.matrix) %in% fireaffected,], diag = TRUE, upper = TRUE)

ecol.ARG.rareREL.FA.d <- vegdist(ecol.ARG.rareREL.matrix[rownames(ecol.ARG.rareREL.matrix) %in% fireaffected,], diag = TRUE, upper = TRUE)

#mantel tests on recovered only
mantel(ecol.rplB.rareREL.FA.d,ecol.ARG.rareREL.FA.d, method = "spear")

#mantel w/ spatial distances
space <- read.table(paste(wd, "/data/spatialdistancematrix.txt", sep = ""), 
                    header=TRUE, row.names=1)

#make spacial matrix names match other data
names(space) <- gsub("C", "cen", names(space))
rownames(space) <- gsub("C", "cen", rownames(space))

#remove rows and columns that involve sites 
#that don't have metagenomes
space <- space[names(space) %in% otu_table$Site, colnames(space) %in% otu_table$Site]

#order spatial distance matrix by other matrices
space.ordered <- space[match(rownames(ecol.rplB.rareREL.matrix), rownames(space)), ]

#make space into a distance matrix
space.d=as.dist(space.ordered, diag = TRUE, upper = TRUE)

#mantel tests v. space
mantel(ecol.rplB.rareREL.d,space.d, method = "spear")
mantel(ecol.ARG.rareREL.d,space.d, method = "spear")

############################
#test rplB and 16S datasets#
############################
#read in 16S OTU table
otu_16s <- fread(paste(wd, "/data/MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", sep = ""), sep = "\t", select = c("OTU_ID",	"C10", "C03", "C14", "C01", "C12", "C06", "C17", "C07", "C04", "C15", "C05", "C16"), colClasses = "numeric")

#make spacial matrix names match other data
names(otu_16s) <- gsub("C", "cen", names(otu_16s))

#transpose data and add rownames
otu_16s_prep <- t(otu_16s)
colnames(otu_16s_prep) <- otu_16s_prep[1,]
otu_16s_prep <- data.frame(otu_16s_prep[-1,])

#make 16S matrix a numeric
otu_16s_prep_num <- data.frame(lapply(otu_16s_prep, function(x) as.numeric(as.character(x))))
rownames(otu_16s_prep_num) <- rownames(otu_16s_prep)

#make 16s data a phyloseq otu table
phyloseq.16s <- otu_table(otu_16s_prep_num, taxa_are_rows = FALSE)

#relativize data
phyloseq.16s.REL <-  transform_sample_counts(phyloseq.16s, function(x) x/sum(x))

#extract 16s matrix 
phyloseq.16s.REL.matrix = as(otu_table(phyloseq.16s.REL), "matrix")

#calculate 16S distance matrix
phyloseq.16s.REL.d <- vegdist(phyloseq.16s.REL.matrix, diag = TRUE, upper = TRUE)

#test 16S v rplB
mantel(ecol.rplB.rareREL.d,phyloseq.16s.REL.d, method = "spear")

#######################
#VENN DIAGRAM ANALYSIS#
#######################
library(VennDiagram)
#prepare data (separate column for soil type)
meta.list <- c("AirTemperature_C","DateFire_Elick2011","OrganicMatter_0360","OrganicMatter_0500","NH4N.ppm_00NA","SulfateSulfur_0ppm","NO3N_0ppm","pH_00NA","K_0ppm","Ca_0ppm","Mg_0ppm","Fe_0ppm","Iron_0ppm", "DateSince_Fire", 'SoilTemperature_to10cm', "As_0ppm")
classification.aes <- aesthetics.blast %>%
  unique(by = "OTU") %>%
  subset(!OTU %in% meta.list) %>%
  subset(Gene !="rplB") %>%
  mutate(PA = (normAbund>0)*1) %>%
  dcast(value.var = "PA", OTU~Classification, fun.aggregate = sum)


grid.newpage()
draw.pairwise.venn(area1 = nrow(subset(classification.aes, FireAffected > 0)), area2 = nrow(subset(classification.aes, Recovered > 0)), cross.area = nrow(subset(classification.aes, Recovered > 0 & FireAffected > 0)), category = c("Fire Affected", "Recovered"))

#calculate total genes assembled (with countable abund)
total.assembled <- (nrow(subset(classification.aes, FireAffected > 0)) + nrow(subset(classification.aes, Recovered > 0)) + nrow(subset(classification.aes, Reference > 0)))

######################################
#TEST OCCURRENCE V ABUNDANCE PATTERNS#
######################################

no <- c("00NA", "0ppm", "0500", "to10cm", "Fire")
abund.occur <- aesthetics.blast %>%
  subset(!number %in% no) %>%
  mutate(PA = (normAbund>0)*1) %>%
  group_by(phylum,class,beginning,OTU) %>%
  summarise(Mean.Perc.Abund = mean(normAbund*100), Site.Occur = sum(PA)) %>%
  subset(Mean.Perc.Abund >0)

#plot abundance v. occurrence
(abund.occur.plot <- ggplot(subset(abund.occur, beginning != "rplB"), aes(x = Site.Occur, y = Mean.Perc.Abund, label = OTU)) +
  geom_jitter(aes(color = beginning), height = 0, width = 0.25, size = 3.5, shape = 1) +
  xlim(0.7,4.3) +
  ylab("Normalized abundance (% rplB)") +
    xlab("Occurrence") +
    scale_color_discrete(name = "Gene") +
      theme_bw(base_size = 14))

#saveplot
ggsave(abund.occur.plot, filename = paste(wd, "/figures/abund.occur.eps", sep = ""), width = 6, height = 4.5)
