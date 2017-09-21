library(vegan)
library(psych)
library(tidyverse)
library(qgraph)
library(phyloseq)
library(reshape2)
library(broom)

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
filenames <- list.files(pattern="*_rformat_dist_0.1.txt")

#move back up directories
setwd("../..")

#make dataframes of all OTU tables
for(i in filenames){
  filepath <- file.path(paste(wd, "/data", sep = ""),paste(i,sep=""))
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
mixed.1 <- c("arxA_12", "arxA_04", "arrA_1", "arrA_2", "arrA_3", "arrA_5", "arrA_6", "arrA_8", "arrA_9", "aioA_18", "aioA_30", "aioA_36", "aioA_44", "aioA_52", "aioA_53", "aioA_57")

#remove OTUs that are gene matches
otu_table <- otu_table[,!names(otu_table) %in% mixed.1]

#rename arxA column that is actually arrA (with no match)
otu_table <- otu_table %>%
  rename(arrA_10 = arxA_10)

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
otu_table_norm_annotated <- data.matrix(otu_table_norm_annotated)
otu_table_norm_annotated.t <- t(otu_table_norm_annotated)

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

#list OTUs present in less than 2 samples
abund <- otu_table_normPA[which(rowSums(otu_table_normPA) > 3),]

#remove OTUs with presence in less than 4 sites
otu_table_norm.slim.t <- t(otu_table_norm_annotated.t[which(rownames(otu_table_norm_annotated.t) %in% rownames(abund)),])

#perform network analysis without rplB!
#remove column based on pattern (rplB)
otu_table_norm.slim.t.genes <- otu_table_norm.slim.t[, -grep("rplB", colnames(otu_table_norm.slim.t))]

#perform correlations
corr.genes <- corr.test(otu_table_norm.slim.t.genes, 
                        method = "spearman", adjust = "fdr", alpha = 0.01)

## prepare network graphics
#read in gene classification data
gene <- read_delim(paste(wd, "/data/gene_classification.txt",  sep=""), 
                   delim = "\t", col_names = TRUE)


#write function to prepare correlation aesthetics 
networkAesthetics <- function(network){
  r.slim <- data.frame(network$r)
  r.slim$gene <- rownames(r.slim)
  r.slim$gene <- gsub("arsC_", "arsC", r.slim$gene)
  r.slim <- r.slim %>%
    separate(gene, c("Gene", "Number"), by = "_", remove = FALSE) %>%
    left_join(gene, by = "Gene") %>%
    select(Gene, Number, Group, gene, gene.color) %>%
    rename(OTU = gene) %>%
    mutate(Shape = "circle")
  r.slim$Group[r.slim$Gene == "rplB"] <- "Organism"
  r.slim$Shape[r.slim$Group == "Organism"] <- "square"
  r.slim$Group[is.na(r.slim$Group)] <- "Metadata"
  r.slim$Shape[r.slim$Group == "Metadata"] <- "diamond"
  r.slim$Shape[r.slim$Group == "AntibioticResistance"] <- "triangle"
  r.final.slim <- r.slim %>%
    separate(OTU, into = c("beginning", "number"), sep = "_")
  r.final.slim$number <- sprintf("%04s", r.final.slim$number)
  r.final.slim$OTU <- paste(r.final.slim$beginning, r.final.slim$number, sep = "_")
  print(r.final.slim)
}

#make network aesthetics 
r.full <- networkAesthetics(corr.genes)

#read in taxize / blast information
blast.final <- read_delim(paste(wd, "/output/taxize_0.1_results_FINAL.txt", sep = ""), col_names = TRUE, delim = " ")

#join aesthetic (shape, taxonomy) information
aesthetics <- r.full %>%
  left_join(blast.final, by = c("OTU", "Gene")) %>%
  unique()

#read in colors for phyla
phylum.colors <- read_delim(paste(wd, "/../networks/data/phylum_colors.txt", sep = ""), delim = "\t", col_names = c("phy.color", "phylum"))

#add phylum colors to aesthetics
aesthetics <- aesthetics %>%
  left_join(phylum.colors, by = "phylum")
aesthetics$phy.color[is.na(aesthetics$gene.color)] <- "#ffffff"

#examine network for initial correlations
clust.network.test <- qgraph(corr.genes$r, minimum = "sig", sampleSize=13, 
                             details = TRUE, layout = "spring",
                             graph = "cor",label.cex = 0.5,
                             alpha = 0.01, graph = "fdr", labels = aesthetics$Gene,  label.scale.equal = TRUE, label.scale = FALSE,shape = aesthetics$Shape, node.resolution = 500,  negDashed = TRUE, curve = 0.2, posCol = "#808080",curveAll = TRUE, overlay = FALSE,  color = aesthetics$gene.color,vsize = 3, GLratio = 6, legend.cex = 0.35, overlay = TRUE)

########################################
#CLUSTER-BASED NETWORK ANALYSIS (+rplB)#
########################################


#List contig clusters in two groups
#from the above network (one neg and one pos
#correlated with temperature)
tneg <- c("ClassB_001", "acr3_036", "acr3_053", "acr3_144", "acr3_002", "acr3_112", "acr3_045", "arsM_034", "arsM_118", "arsM_296", "arsM_109", "arsM_389", "arsC_glut_103", "arsC_glut_177", "arsC_glut_126", "tolC_08", "dfra12_076", "dfra12_038", "dfra12_074", "dfra12_058", "dfra12_046", "dfra12_073", "vanX_03", "arsA_011")
tpos <- c("vanX_17", "arsA_054", "arsA_056", "arsA_012", "arsA_063", "arsA_033", "arsA_005", "arsA_022", "vanA_10", "arsC_glut_227", "arsC_glut_182", "arsC_glut_112", "arsC_glut_054", "arsC_glut_071", "arsC_glut_044", "arsC_glut_043", "arsM_188", "arsM_198", "arsM_060", "arsM_201", "arsM_023", "arsM_189", "arsM_047", "arsM_465", "arsM_159", "arsM_009", "arsM_387", "arsM_275", "arsM_182", "arsM_101", "arsM_262", "arsM_372", "arsM_073", "arsM_072", "arsM_120", "sul2_01", "dfra12_081", "acr3_078", "acr3_031", "acr3_065", "acr3_160", "acr3_047", "acr3_094", "acr3_127", "acr3_152", "acr3_048", "acr3_143", "acr3_090", "acr3_120", "acr3_101", "acr3_018", "acr3_106", "ClassA_01", "ClassA_06")

#extract rplB sequences and cluster seqs
tneg_otu_table <- otu_table_norm.slim.t[,colnames(otu_table_norm.slim.t) %in% tneg]
tpos_otu_table <- otu_table_norm.slim.t[,colnames(otu_table_norm.slim.t) %in% tpos]
rplB_otu_table <- otu_table_norm.slim.t[, grep("rplB", colnames(otu_table_norm.slim.t))]

#join rplB data with cluster seqs
tneg_otu_table.lrg <- cbind(tneg_otu_table, rplB_otu_table)  
tpos_otu_table.lrg <- cbind(tpos_otu_table, rplB_otu_table)  

#perform correlations
tneg.net <- corr.test(tneg_otu_table.lrg, 
                      method = "spearman", adjust = "fdr", alpha = 0.01)
tpos.net <- corr.test(tpos_otu_table.lrg, 
                      method = "spearman", adjust = "fdr", alpha = 0.01)

#extract rplB significantly correlated with clusters
tneg.net.rplB <- tneg.net$p %>%
  melt() %>%
  rename(Contig1 = Var1, Contig2 = Var2, p = value) %>%
  subset(p < 0.01)

tneg.net.rplB <- tneg.net.rplB[!grepl("rplB", tneg.net.rplB$Contig2),]
tneg.net.rplB <- tneg.net.rplB[grepl("rplB", tneg.net.rplB$Contig1),]

tpos.net.rplB <- tpos.net$p %>%
  melt() %>%
  rename(Contig1 = Var1, Contig2 = Var2, p = value) %>%
  subset(p < 0.01)

tpos.net.rplB <- tpos.net.rplB[!grepl("rplB", tpos.net.rplB$Contig2),]
tpos.net.rplB <- tpos.net.rplB[grepl("rplB", tpos.net.rplB$Contig1),]

#trim contig tables based on co-occurring rplB
tpos.final <- cbind(tpos_otu_table, rplB_otu_table[,colnames(rplB_otu_table) %in% tpos.net.rplB$Contig1])
tneg.final <- cbind(tneg_otu_table, rplB_otu_table[,colnames(rplB_otu_table) %in% tneg.net.rplB$Contig1])

##perform SLIM correlations
tneg.net.final <- corr.test(tneg.final, 
                            method = "spearman", adjust = "fdr", alpha = 0.01)
tpos.net.final <- corr.test(tpos.final, 
                            method = "spearman", adjust = "fdr", alpha = 0.01)

#make network aesthetics 
tneg.r <- networkAesthetics(tneg.net.final)
tpos.r <- networkAesthetics(tpos.net.final)

#join aesthetic (shape, taxonomy) information
aesthetics.tneg <- tneg.r %>%
  left_join(blast.final, by = c("OTU", "Gene")) %>%
  unique()

aesthetics.tpos <- tpos.r %>%
  left_join(blast.final, by = c("OTU", "Gene")) %>%
  unique()

#add phylum colors to aesthetics
aesthetics.tneg <- aesthetics.tneg %>%
  left_join(phylum.colors, by = "phylum")
aesthetics.tneg$phy.color[is.na(aesthetics.tneg$phy.color)] <- "#ffffff"

aesthetics.tpos <- aesthetics.tpos %>%
  left_join(phylum.colors, by = "phylum")
aesthetics.tpos$phy.color[is.na(aesthetics.tpos$phy.color)] <- "#ffffff"

#visualize networks
tneg.network <- qgraph(tneg.net.final$r, minimum = "sig", sampleSize=13, 
                       details = TRUE, layout = "spring",
                       graph = "cor",label.cex = 0.7,
                       alpha = 0.01, graph = "fdr", label.scale.equal = TRUE, label.scale = FALSE, negDashed = TRUE, curve = 0.2, posCol = "#808080",curveAll = TRUE, shape = aesthetics.tneg$Shape, color = aesthetics.tneg$phy.color, labels = aesthetics.tneg$Gene,vsize = 3)

tpos.network <- qgraph(tpos.net.final$r, minimum = "sig", sampleSize=13, 
                       details = TRUE, layout = "spring",
                       graph = "cor",label.cex = 0.7,
                       alpha = 0.01, graph = "fdr", label.scale.equal = TRUE, label.scale = FALSE,negDashed = TRUE, curve = 0.2, posCol = "#808080",curveAll = TRUE, shape = aesthetics.tpos$Shape, color= aesthetics.tpos$phy.color, vsize = 3.2, labels = aesthetics.tpos$Gene)

#now plot abundance/ location of these gene clusters
tpos.abund <- otu_table_norm_annotated.t[rownames(otu_table_norm_annotated.t) %in% colnames(tpos.final),]
tpos.abund.annotated <- tpos.abund %>%
  melt() %>% rename(OTU = Var1, Site = Var2, Abundance = value) %>%
  left_join(meta, by = "Site") %>%
  separate(OTU, into = c("Gene", "Number"))

tpos.abund.annotated$Site <- factor(tpos.abund.annotated$Site, 
                                    tpos.abund.annotated$Site[order(tpos.abund.annotated$SoilTemperature_to10cm)])
(tpos.abund <- ggplot(subset(tpos.abund.annotated, Gene !="rplB"), aes(x = Site, y = Abundance)) +
    geom_bar(stat = "identity", color = "black", aes(fill = Gene)) +
    theme_bw(base_size = 10) +
    scale_fill_manual(values = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#CCCCCC", "#FFFFCC", "#B3E2CD", "#FFFFFF", "#F2F2F2")) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=0.95,vjust=0.9)))

ggsave(tpos.abund, filename = paste(wd, "/figures/abundance_tpos.eps", sep= ""), units = "in", width = 4, height = 3)

(tpos.abund.rplB <- ggplot(subset(tpos.abund.annotated, Gene =="rplB"), aes(x = Site, y = Abundance)) +
    geom_bar(stat = "identity", color = "black", fill = "#E3D7BC") +
    theme_bw(base_size = 10) +
    ylim(0, 0.73) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=0.95,vjust=0.9)))

ggsave(tpos.abund.rplB, filename = paste(wd, "/figures/abundance_tpos_rplB.eps", sep= ""), units = "in", width = 3, height = 3)


#repeat for neg temp cluster
tneg.abund <- otu_table_norm_annotated.t[rownames(otu_table_norm_annotated.t) %in% colnames(tneg.final),]
tneg.abund.annotated <- tneg.abund %>%
  melt() %>% rename(OTU = Var1, Site = Var2, Abundance = value) %>%
  left_join(meta, by = "Site") %>%
  separate(OTU, into = c("Gene", "Number"))

tneg.abund.annotated$Site <- factor(tneg.abund.annotated$Site, 
                                    tneg.abund.annotated$Site[order(tneg.abund.annotated$SoilTemperature_to10cm)])
(tneg.abund <- ggplot(subset(tneg.abund.annotated, Gene !="rplB"), aes(x = Site, y = Abundance)) +
    geom_bar(stat = "identity", color = "black", aes(fill = Gene)) +
    theme_bw(base_size = 10) +
    scale_fill_manual(values = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#FDDAEC", "#F2F2F2")) +
    ylim(0, 0.25) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=0.95,vjust=0.9)))
ggsave(tneg.abund, filename = paste(wd, "/figures/abundance_tneg.eps", sep= ""), units = "in", width = 4, height = 3)

(tneg.abund.rplB <- ggplot(subset(tneg.abund.annotated, Gene =="rplB"), aes(x = Site, y = Abundance)) +
    geom_bar(stat = "identity", color = "black", fill = "#E3D7BC") +
    theme_bw(base_size = 10) +
    ylim(0, 0.25) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=0.95,vjust=0.9)))
ggsave(tneg.abund.rplB, filename = paste(wd, "/figures/abundance_tneg_rplB.eps", sep= ""), units = "in", width = 3, height = 3)


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

#order based on group
gene_abundance_summary$Gene <- factor(gene_abundance_summary$Gene, 
                                      levels = gene_abundance_summary$Gene[order(gene_abundance_summary$Description)])

#plot antibiotic resistance genes
(boxplotARG <- ggplot(subset(gene_abundance_summary, subset = Group == "AntibioticResistance"), aes(x = Classification, 
                                                                                                    y = Total*100)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm)) +
    facet_wrap(~Gene, scales = "free_y", ncol = 3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Gene per rplB (%)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=0.95)))

ggsave(boxplotARG, filename = paste(wd, "/figures/ARGboxplot.eps", sep = ""), height = 8, width = 7, units = "in")

(boxplot.asrg <- ggplot(subset(gene_abundance_summary, subset = Group == "ArsenicResistance"), aes(x = Classification, 
                                                                                                   y = Total*100)) +
    geom_boxplot() +
    geom_jitter(aes(color = SoilTemperature_to10cm)) +
    facet_wrap(~Gene, scales = "free_y", ncol = 3) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    ylab("Gene per rplB (%)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=0.95)))

ggsave(boxplot.asrg, filename = paste(wd, "/figures/AsRGboxplot.eps", sep = ""), height = 5, width = 7, units = "in")

######################
#MANN WHITNEY U TESTS#
######################

#perform test of gene abundance compared to 
#soil history (Classification)
mwu.classification <- subset(gene_abundance_summary, Classification !="Reference") %>% group_by(Gene) %>% do(tidy(wilcox.test(.$Total~.$Classification, paired = FALSE)))

#save table
write.table(mwu.classification, paste(wd, "/output/mannwhitneyu.csv", sep = ""), row.names = FALSE, sep = ",", quote = FALSE)

######################
#CORRELATION ANALYSES#
######################

#decast for abundance check but first make Site a character
gene_abundance_summary$Site <- as.character(gene_abundance_summary$Site)
cast.gene <- data.frame(acast(gene_abundance_summary, Site ~ Gene, value.var = "Total"))

#call na's zeros
cast.gene[is.na(cast.gene)] =0

#remove unnecessary sites in meta data
meta.slim <- meta[meta$Site %in% gene_abundance_summary$Site,]

#correlate data with temp
corr <- print(corr.test(x = cast.gene, y = meta.slim[,c(2,16, 20:30)], method = "spearman", adjust = "fdr"), short = FALSE)
cor <- corr.test(x = cast.gene, y = meta.slim[,c(2,16, 20:30)], method = "spearman", adjust = "fdr")
#arsenic and antibiotic resistance genes are not correlated with temperature!

#save as table for supplemental material
write.table(corr, file = paste(wd, "/output/geochem.correlations.csv", sep = ""), sep = ",", row.names = TRUE, quote = FALSE)

#check correlations between genes
cast.gene <- cast.gene %>% select(-c(rplB, arrA, arxA, tetA, tetW, tetX))
cast.gene <- cast.gene[,-1]
corr.genes <- print(corr.test(cast.gene, method = "spearman", adjust = "fdr"), 
                    short = FALSE)
corr.genes.matrix <- corr.test(cast.gene, method = "spearman", adjust = "fdr")

#save correlations to table
write.table(corr.genes, paste(wd, "/output/gene.correlations.csv", sep = ""), row.names = TRUE, sep = ",", quote = FALSE)

cor.plot(corr.genes.matrix$r,numbers = TRUE, xlas = 2, upper = FALSE, diag = FALSE, stars = TRUE, pval = corr.genes.matrix$p)


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
setwd("../")
wd <- print(getwd())

#split columns and tidy dataset
data <- data %>%
  separate(col = id, into = c("Site", "junk"), sep = 5, remove = TRUE) %>%
  separate(col = junk, into = c("Gene", "junk"), sep = "_45_", remove = TRUE)

#remove awkward _ in Gene column
data$Gene <- gsub("_", "", data$Gene)

#change site from "cen" to "Cen" so it matches metadata
data$Site <- gsub("cen", "Cen", data$Site)

#separage out rplB data (not needed for gene-centric analysis)
rplB <- data[which(data$Gene == "rplB"),]
data <- data[-which(data$Gene == "rplB"),]

#split columns 
rplB <- rplB %>%
  select(Site, Taxon:Fraction.Abundance) %>%
  group_by(Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
rplB$Fraction.Abundance <- as.numeric(rplB$Fraction.Abundance)
rplB$Abundance <- as.numeric(rplB$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.rplB <- rplB %>%
  summarise(Total = sum(Fraction.Abundance), rplB = sum(Abundance))

#decast for abundance check and call na's zero
dcast <- acast(rplB, Taxon ~ Site, value.var = "Fraction.Abundance")
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

#########################################################
#a-DIVERSITY ANALYSIS OF ARGs AND AsRGs clustered at 0.1#
#########################################################

#separate original (non-normalized) contig table into
#specific gene groups (rplB, ARG, AsRG)
ecol.rplB <- otu_table[, grep("rplB", colnames(otu_table))]
ecol.AsRG <- otu_table[, grep("ars|aio|arx|arr|acr", 
                              colnames(otu_table))]
ecol.ARG <- otu_table[,grep("tolC|dfra|Class|tet|van|CEP|AAC|ade|sul", 
                            colnames(otu_table))]

#add site names as row names to each contig table 
rownames(ecol.rplB) <- otu_table[,1]
rownames(ecol.ARG) <- otu_table[,1]
rownames(ecol.AsRG) <- otu_table[,1]

#make list of 13 colors (based on classification)
class.13 <- c("#FFFF00", "darkgreen", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000")

#check sampling depth of each matrix
rarecurve(ecol.rplB, step=1, label = FALSE, col = class.13)
rarecurve(ecol.AsRG, step=1, label = FALSE, col = class.13)
rarecurve(ecol.ARG, step=1, label = FALSE, col = class.13)

#convert contig tables into phyloseq objects
ecol.rplB.phyloseq <- otu_table(ecol.rplB, taxa_are_rows = FALSE)
ecol.AsRG.phyloseq <- otu_table(ecol.AsRG, taxa_are_rows = FALSE)
ecol.ARG.phyloseq <- otu_table(ecol.ARG, taxa_are_rows = FALSE)

#rarefy to even sampling depth 
ecol.rplB.rare <- rarefy_even_depth(ecol.rplB.phyloseq, rngseed = TRUE)
rarecurve(ecol.rplB.rare, step=1, label = FALSE, col = class.13)

ecol.ARG.rare <- rarefy_even_depth(ecol.ARG.phyloseq, rngseed = TRUE)
rarecurve(ecol.ARG.rare, step=1, label = FALSE, col = class.13)

ecol.AsRG.rare <- rarefy_even_depth(ecol.AsRG.phyloseq, rngseed = TRUE)
rarecurve(ecol.AsRG.rare, step=1, label = FALSE, col = class.13)

#calculate evenness
plieou.rplB <- data.frame(group = "rplB", Site = rownames(ecol.rplB.rare), plieou = vegan::diversity(ecol.rplB.rare, index = "shannon")/log(specnumber(ecol.rplB.rare)))

plieou.ARG <- data.frame(group = "ARG", Site = rownames(ecol.ARG.rare), plieou = vegan::diversity(ecol.ARG.rare, index = "shannon")/log(specnumber(ecol.ARG.rare)))

plieou.AsRG <- data.frame(group = "AsRG", Site = rownames(ecol.AsRG.rare), plieou = vegan::diversity(ecol.AsRG.rare, index = "shannon")/log(specnumber(ecol.AsRG.rare)))

#join all evenness information and add metadata
plieou.full <- rbind(plieou.ARG, plieou.AsRG, plieou.rplB)
plieou.full$Site <- gsub("cen", "Cen", plieou.full$Site)
plieou.full <- left_join(plieou.full, meta, by = "Site")

#plot evenness
(plieou.plot <- ggplot(plieou.full, aes(x = Classification, y = plieou)) +
    geom_boxplot() +
    ylab(label = "Evenness") +
    theme_bw() +
    facet_wrap(~group) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

#save evennes plots
ggsave(plieou.plot, filename = paste(wd, "/figures/evenness.eps", sep = ""))

#test differences in evenness
#mann whitney u test for significance
#perform test of gene abundance compared to 
#soil history (Classification)
evenness.mwu.classification <- subset(plieou.full, Classification !="Reference") %>% group_by(group) %>% do(tidy(wilcox.test(abs(.$plieou)~.$Classification, paired = FALSE)))

#save table
write.table(evenness.mwu.classification, paste(wd, "/output/evenness.1.mannwhitneyu.csv", sep = ""), row.names = FALSE, sep = ",", quote = FALSE)

#make metadata a phyloseq class object
meta$Site <- gsub("Cen", "cen", meta$Site)
rownames(meta) <- meta[,1]
meta.phylo <- meta[,-1]
meta.phylo <- sample_data(meta.phylo)

#read in trees and make phyloseq object
#library(ape)
#rplb.tree <- read.tree(paste(wd, "/data/rplB_0.1_FastTree.nwk", sep = ""))
#rplb.tree <- phy_tree(rplb.tree)

##make biom for phyloseq
ecol.rplB.rare <- merge_phyloseq(ecol.rplB.rare, meta.phylo)
ecol.ARG.rare <- merge_phyloseq(ecol.ARG.rare, meta.phylo)
ecol.AsRG.rare <- merge_phyloseq(ecol.AsRG.rare, meta.phylo)

#plot & save RICHNESS
(richness.ecol.rplB.rare <- plot_richness(ecol.rplB.rare, x = "Classification", measures = "Observed") +
    geom_boxplot() +
    theme_bw() +
    ylim(30,250) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

(richness.ecol.ARG.rare <- plot_richness(ecol.ARG.rare, x = "Classification", measures = "Observed") +
    geom_boxplot() +
    theme_bw() +
    ylim(30,250) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

(richness.ecol.AsRG.rare <- plot_richness(ecol.AsRG.rare, x = "Classification", measures = "Observed") +
    geom_boxplot() +
    theme_bw() +
    ylim(30,250) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

#extract richness and perform statistical tests
richness.ecol.rplB.rare.data <- estimate_richness(ecol.rplB.rare, split = TRUE, measures = "Observed")
richness.ecol.rplB.rare.data <- richness.ecol.rplB.rare.data %>%
  mutate(Site = rownames(richness.ecol.rplB.rare.data), GeneGroup = "rplB")

richness.ecol.ARG.rare.data <- estimate_richness(ecol.ARG.rare, split = TRUE, measures = "Observed")
richness.ecol.ARG.rare.data <- richness.ecol.ARG.rare.data %>%
  mutate(Site = rownames(richness.ecol.ARG.rare.data), GeneGroup = "ARG")

richness.ecol.AsRG.rare.data <- estimate_richness(ecol.AsRG.rare, split = TRUE, measures = "Observed")
richness.ecol.AsRG.rare.data <- richness.ecol.AsRG.rare.data %>%
  mutate(Site = rownames(richness.ecol.AsRG.rare.data), GeneGroup = "AsRG")

#join together richness results and add metadata 
richness.data <- rbind(richness.ecol.rplB.rare.data, richness.ecol.ARG.rare.data, richness.ecol.AsRG.rare.data)
richness.data <- richness.data %>%
  left_join(meta, by = "Site")

#mann whitney u test for significance
#perform test of gene abundance compared to 
#soil history (Classification)
richness.mwu.classification <- subset(richness.data, Classification !="Reference") %>% group_by(GeneGroup) %>% do(tidy(wilcox.test(abs(.$Observed)~.$Classification, paired = FALSE)))

#save table
write.table(richness.mwu.classification, paste(wd, "/output/richness.1.mannwhitneyu.csv", sep = ""), row.names = FALSE, sep = ",", quote = FALSE)

#relativize rarefied datasets before beta-diversity
ecol.rplB.rareREL <-  transform_sample_counts(ecol.rplB.rare, function(x) x/sum(x))
ecol.ARG.rareREL <-  transform_sample_counts(ecol.ARG.rare, function(x) x/sum(x))
ecol.AsRG.rareREL <- transform_sample_counts(ecol.AsRG.rare, function(x) x/sum(x))

#########################################################
#B-DIVERSITY ANALYSIS OF ARGs AND AsRGs clustered at 0.1#
#########################################################

#plot Bray Curtis ordination for rplB
ord.rplB.bray <- ordinate(ecol.rplB.rareREL, method="PCoA", distance="bray")
(bc.ord.rplB=plot_ordination(ecol.rplB.rareREL, ord.rplB.bray, color="SoilTemperature_to10cm",
                             title="Bray Curtis", shape = "Classification") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    geom_point(size=5) +
    theme_light(base_size = 12))

#plot Bray Curtis ordination for ARGs
ord.ARG.bray <- ordinate(ecol.ARG.rareREL, method="PCoA", distance="bray")
(bc.ord.ARG=plot_ordination(ecol.ARG.rareREL, ord.ARG.bray, color="SoilTemperature_to10cm",
                            title="Bray Curtis", shape = "Classification") +
    geom_point(size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#plot Bray Curtis ordination for AsRGs
ord.AsRG.bray <- ordinate(ecol.AsRG.rareREL, method="PCoA", distance="bray")
(bc.ord.AsRG=plot_ordination(ecol.AsRG.rareREL, ord.AsRG.bray, color="SoilTemperature_to10cm",
                             title="Bray Curtis", shape = "Classification") +
    geom_point(size=5) +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    theme_light(base_size = 12))

#save each bray curtis
ggsave(bc.ord.rplB, filename = paste(wd, "/figures/rplB.braycurtis.eps", sep = ""), width = 5, height =3, units = "in")
ggsave(bc.ord.ARG, filename = paste(wd, "/figures/ARG.braycurtis.eps", sep = ""), width = 5, height = 3, units = "in")
ggsave(bc.ord.AsRG, filename = paste(wd, "/figures/AsRG.braycurtis.eps", sep = ""), width = 5, height = 3, units = "in")

#extract OTU table from phyloseq
ecol.rplB.rareREL.matrix = as(otu_table(ecol.rplB.rareREL), "matrix")
ecol.AsRG.rareREL.matrix = as(otu_table(ecol.AsRG.rareREL), "matrix")
ecol.ARG.rareREL.matrix = as(otu_table(ecol.ARG.rareREL), "matrix")

#calculate distance matrix
ecol.rplB.rareREL.d <- vegdist(ecol.rplB.rareREL.matrix, diag = TRUE, upper = TRUE)
ecol.ARG.rareREL.d <- vegdist(ecol.ARG.rareREL.matrix, diag = TRUE, upper = TRUE)
ecol.AsRG.rareREL.d <- vegdist(ecol.AsRG.rareREL.matrix, diag = TRUE, upper = TRUE)

#mantel tests
mantel(ecol.rplB.rareREL.d,ecol.ARG.rareREL.d, method = "spear")
mantel(ecol.rplB.rareREL.d,ecol.AsRG.rareREL.d, method = "spear")
mantel(ecol.AsRG.rareREL.d,ecol.ARG.rareREL.d, method = "spear")

#try mantel tests of recovered sites only
recovered <- c("cen01", "cen03", "cen04", "cen05", "cen07")

ecol.rplB.rareREL.Rec.d <- vegdist(ecol.rplB.rareREL.matrix[rownames(ecol.rplB.rareREL.matrix) %in% recovered,], diag = TRUE, upper = TRUE)

ecol.ARG.rareREL.Rec.d <- vegdist(ecol.ARG.rareREL.matrix[rownames(ecol.ARG.rareREL.matrix) %in% recovered,], diag = TRUE, upper = TRUE)

ecol.AsRG.rareREL.Rec.d <- vegdist(ecol.AsRG.rareREL.matrix[rownames(ecol.AsRG.rareREL.matrix) %in% recovered,], diag = TRUE, upper = TRUE)

#mantel tests on recovered only
mantel(ecol.rplB.rareREL.Rec.d,ecol.ARG.rareREL.Rec.d, method = "spear")
mantel(ecol.rplB.rareREL.Rec.d,ecol.AsRG.rareREL.Rec.d, method = "spear")
mantel(ecol.AsRG.rareREL.Rec.d,ecol.ARG.rareREL.Rec.d, method = "spear")

#try mantel tests on fire affected sites only
fireaffected <- c("cen06", "cen10", "cen12", "cen13", "cen14", "cen15", "cen16")

ecol.rplB.rareREL.FA.d <- vegdist(ecol.rplB.rareREL.matrix[rownames(ecol.rplB.rareREL.matrix) %in% fireaffected,], diag = TRUE, upper = TRUE)

ecol.ARG.rareREL.FA.d <- vegdist(ecol.ARG.rareREL.matrix[rownames(ecol.ARG.rareREL.matrix) %in% fireaffected,], diag = TRUE, upper = TRUE)

ecol.AsRG.rareREL.FA.d <- vegdist(ecol.AsRG.rareREL.matrix[rownames(ecol.AsRG.rareREL.matrix) %in% fireaffected,], diag = TRUE, upper = TRUE)

#mantel tests on recovered only
mantel(ecol.rplB.rareREL.FA.d,ecol.ARG.rareREL.FA.d, method = "spear")
mantel(ecol.rplB.rareREL.FA.d,ecol.AsRG.rareREL.FA.d, method = "spear")
mantel(ecol.AsRG.rareREL.FA.d,ecol.ARG.rareREL.FA.d, method = "spear")

#mantel w/ spatial distances
space <- read.table(paste(wd, "/data/spatialdistancematrix.txt", sep = ""), 
                    header=TRUE, row.names=1)

#make spacial matrix names match other data
names(space) <- gsub("C", "cen", names(space))
rownames(space) <- gsub("C", "cen", rownames(space))

#remove rows and columns that involve sites 
#that don't have metagenomes
space <- space[names(space) %in% otu_table$Site, colnames(space) %in% otu_table$Site]

#make space into a distance matrix
space.d=as.dist(space, diag = TRUE, upper = TRUE)

#mantel tests v. space
mantel(ecol.rplB.rareREL.d,space.d, method = "spear")
mantel(ecol.ARG.rareREL.d,space.d, method = "spear")
mantel(ecol.AsRG.rareREL.d,space.d, method = "spear")

####################
#TAXANOMIC ANALYSIS#
####################

#add taxanomic informaiton (BLAST) to gene_abundance files
#first need to add leading zeros to gene abundance OTUs
gene_abundance$OTU <- gsub("arsC_", "arsC", gene_abundance$OTU)
gene_abundance_blast <- gene_abundance %>%
  separate(OTU, into = c("prefix", "number"))
gene_abundance_blast$number <- sprintf("%04s", gene_abundance_blast$number)
gene_abundance_blast$number <- paste(gene_abundance_blast$prefix, gene_abundance_blast$number, sep = "_")

#join with BLAST data
gene_abundance_blast_final <- gene_abundance_blast %>%
  rename(OTU = number) %>%
  left_join(blast.final, by = "OTU") %>%
  group_by(Classification, Site, Group, phylum) %>%
  summarise(Abund = sum(RelativeAbundance)) %>%
  ungroup() %>%
  group_by(Classification, Group, phylum) %>%
  summarise(meanAbund = mean(Abund))

#set rplB as Community membership (rplB)
gene_abundance_blast_final$Group[is.na(gene_abundance_blast_final$Group)] <- "rplB"

#select top 10 phyla (from previous rplB analysis)
top.phy <- c("Acidobacteria", "Actinobacteria","Bacteroidetes",  "Chloroflexi", "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Verrucomicrobia")

#plot based on taxonomy (separate rplB)
ggplot(subset(gene_abundance_blast_final, phylum %in% top.phy), aes(x = phylum, y = meanAbund, fill = Group, group = Group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  facet_wrap(~Classification) +
  theme(axis.text.x = element_text(angle = 90, size = 8, 
                                   hjust=0.95,vjust=0.2))

gene_abundance_cast <- dcast(gene_abundance_blast_final, Classification+phylum ~ Group, fun.aggregate = sum)
la <- subset(gene_abundance_cast, phylum %in% top.phy)

gene_abundance_cast %>% group_by(phylum) %>% do(tidy(cor.test(.$rplB, .$AntibioticResistance)))

ggplot(subset(gene_abundance_cast, phylum %in% top.phy), aes(x = rplB, y = ArsenicResistance, color = phylum)) +
  geom_point()
