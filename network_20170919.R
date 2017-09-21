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
write.table(rplB_summary_save, file = paste(wd, "/output/rplB.summary.scg.txt", sep = ""), row.names = FALSE)


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
abund <- otu_table_normPA[which(rowSums(otu_table_normPA) > 9),]

#remove OTUs with presence in less than 4 sites
otu_table_norm.slim.t <- t(otu_table_norm_annotated.t[which(rownames(otu_table_norm_annotated.t) %in% rownames(abund)),])

ecol.AsRG <- otu_table_norm.slim.t[, grep("ars|aio|arx|arr|acr", 
                              colnames(otu_table_norm.slim.t))]
ecol.ARG <- otu_table_norm.slim.t[,grep("tolC|dfra|Class|tet|van|CEP|AAC|ade|sul", 
                            colnames(otu_table_norm.slim.t))]

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

#examine network for initial correlations
clust.network.test <- qgraph(corr.genes$r, minimum = "sig", sampleSize=13, 
                        details = TRUE, layout = "spring",
                        graph = "cor",label.cex = 0.5,
                        alpha = 0.01, graph = "fdr", labels = aesthetics$OTU,  label.scale.equal = TRUE, label.scale = FALSE,shape = aesthetics$Shape, node.resolution = 500,  negDashed = TRUE, curve = 0.2, posCol = "#808080",curveAll = TRUE, overlay = FALSE,  color = aesthetics$gene.color,vsize = 4, GLratio = 6, legend.cex = 0.35, overlay = TRUE)

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

