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
#setwd("/Users/dunivint/Documents/GitHubRepos/ARG-AsRG_co-occurrence_Centralia/diversity_analysis/")
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
filenames <- list.files(pattern="*_rformat_dist_0.01.txt")

#move back up directories
setwd("../")

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
abund <- otu_table_normPA[which(rowSums(otu_table_normPA) > 2),]

#remove OTUs with presence in less than 4 sites
otu_table_norm.slim.t <- t(otu_table_norm_annotated.t[which(rownames(otu_table_norm_annotated.t) %in% rownames(abund)),])

#perform network analysis without rplB!
#remove column based on pattern (rplB)
otu_table_norm.slim.t.genes <- otu_table_norm.slim.t[, -grep("rplB", colnames(otu_table_norm.slim.t))]

#perform correlations
corr.genes <- corr.test(otu_table_norm.slim.t, 
                        method = "spearman", adjust = "fdr", alpha = 0.05)

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
blast.final <- read_delim(paste(wd, "/output/taxize_0.01_results_FINAL.txt", sep = ""), col_names = TRUE, delim = " ")

#join aesthetic (shape, taxonomy) information
aesthetics <- r.full %>%
  left_join(blast.final, by = c("OTU", "Gene")) %>%
  unique()

#read in colors for phyla
phylum.colors <- read_delim(paste(wd, "/data/phylum_colors.txt", sep = ""), delim = "\t", col_names = c("phy.color", "phylum"))

#add phylum colors to aesthetics
aesthetics <- aesthetics %>%
  left_join(phylum.colors, by = "phylum")
aesthetics$phy.color[is.na(aesthetics$phylum)] <- "#ffffff"

#examine network for initial correlations
clust.network.test <- qgraph(corr.genes$r, minimum = "sig", sampleSize=13, 
                             details = TRUE, layout = "spring",
                             graph = "cor",label.cex = 0.6,
                             alpha = 0.05, graph = "fdr", labels = aesthetics$Gene,  label.scale.equal = TRUE, label.scale = FALSE,shape = aesthetics$Shape, node.resolution = 500,  negDashed = TRUE, curve = 0.2, posCol = "#808080",curveAll = TRUE, color = aesthetics$phy.color,overlay = FALSE,  vsize = 3, overlay = TRUE)

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
corr <- print(corr.test(x = cast.gene, y = meta.slim[,c(2,16, 20:30)], method = "spearman", adjust = "fdr"), short = FALSE)
cor <- corr.test(x = cast.gene, y = meta.slim[,c(2,16, 20:30)], method = "spearman", adjust = "fdr")
#arsenic and antibiotic resistance genes are not correlated with temperature!

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

sig <- c("ClassA", "ClassB", "dfra12","sul2", "tolC")
#plot antibiotic resistance genes
(sigTempGene <- ggplot(subset(gene_abundance_summary, subset = Gene %in% sig), aes(x = log(SoilTemperature_to10cm), 
                                                                                                    y = Total)) +
    geom_point(aes(shape = Classification), size = 2) +
    #geom_jitter(aes(color = SoilTemperature_to10cm)) +
    facet_wrap(~Gene, scales = "free_y", ncol = 2) +
    ylab("rplB-normalized abundance") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=0.95)))

ggsave(sigTempGene, filename = paste(wd, "/figures/sig_temp_gene.png", sep = ""), height = 8, width = 7, units = "in")

(nonTempGene <- ggplot(subset(gene_abundance_summary, subset = !Gene %in% sig), aes(x = log(SoilTemperature_to10cm), 
                                                                                   y = Total)) +
    geom_point(aes(shape = Classification), size = 2) +
    #geom_jitter(aes(color = SoilTemperature_to10cm)) +
    facet_wrap(~Gene, scales = "free_y", ncol = 3) +
    ylab("rplB-normalized abundance") +
    theme_bw(base_size = 12))

ggsave(nonTempGene, filename = paste(wd, "/figures/nonTempGene.png", sep = ""), height = 8, width = 7, units = "in")

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

#remove Cen13
rplB <- rplB[!rplB$Site == "Cen13",]

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
#a-DIVERSITY ANALYSIS OF ARGs AND AsRGs clustered at 0.01#
#########################################################

#separate original (non-normalized) contig table into
#specific gene groups (rplB, ARG, AsRG)
ecol.rplB <- otu_table[, grep("rplB", colnames(otu_table))]
ecol.ARG <- otu_table[,grep("tolC|dfra|Class|tet|van|CEP|AAC|ade|sul", 
                            colnames(otu_table))]

#add site names as row names to each contig table 
rownames(ecol.rplB) <- otu_table[,1]
rownames(ecol.ARG) <- otu_table[,1]

#make list of 13 colors (based on classification)
class.13 <- c("#FFFF00", "darkgreen", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FFFF00", "#FF0000")

#check sampling depth of each matrix
rarecurve(ecol.rplB, step=1, label = FALSE, col = class.13)
rarecurve(ecol.ARG, step=1, label = FALSE, col = class.13)

#convert contig tables into phyloseq objects
ecol.rplB.phyloseq <- otu_table(ecol.rplB, taxa_are_rows = FALSE)
ecol.ARG.phyloseq <- otu_table(ecol.ARG, taxa_are_rows = FALSE)

#rarefy to even sampling depth 
ecol.rplB.rare <- rarefy_even_depth(ecol.rplB.phyloseq, rngseed = TRUE)
rarecurve(ecol.rplB.rare, step=1, label = FALSE, col = class.13)

ecol.ARG.rare <- rarefy_even_depth(ecol.ARG.phyloseq, rngseed = TRUE)
rarecurve(ecol.ARG.rare, step=1, label = FALSE, col = class.13)

#calculate evenness
plieou.rplB <- data.frame(group = "rplB", Site = rownames(ecol.rplB.rare), plieou = vegan::diversity(ecol.rplB.rare, index = "shannon")/log(specnumber(ecol.rplB.rare)))

plieou.ARG <- data.frame(group = "ARG", Site = rownames(ecol.ARG.rare), plieou = vegan::diversity(ecol.ARG.rare, index = "shannon")/log(specnumber(ecol.ARG.rare)))


#join all evenness information and add metadata
plieou.full <- rbind(plieou.ARG, plieou.rplB)
plieou.full$Site <- gsub("cen", "Cen", plieou.full$Site)
plieou.full <- left_join(plieou.full, meta, by = "Site")

#plot evenness
(plieou.plot <- ggplot(plieou.full, aes(x = SoilTemperature_to10cm, y = plieou)) +
    geom_point(aes(shape = Classification), size = 2) +
    ylab(label = "Evenness") +
    theme_bw(base_size = 12) +
    facet_wrap(~group))

#save evennes plots
ggsave(plieou.plot, filename = paste(wd, "/figures/evenness.eps", sep = ""))

#correlate evenness with Temp
spearman_evenness <- plieou.full %>% group_by(group) %>% do(tidy(cor.test(.$plieou, .$SoilTemperature_to10cm, method = "spearman")))

#save as table for supplemental material
write.table(spearman_evenness, file = paste(wd, "/output/evenness.Temp.correlations.csv", sep = ""), sep = ",", row.names = FALSE, quote = FALSE)

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

##make biom for phyloseq
ecol.rplB.rare <- merge_phyloseq(ecol.rplB.rare, meta.phylo)
ecol.ARG.rare <- merge_phyloseq(ecol.ARG.rare, meta.phylo)

#plot & save RICHNESS
(richness.ecol.rplB.rare <- plot_richness(ecol.rplB.rare, x = "SoilTemperature_to10cm", measures = "Observed") +
    geom_point(aes(shape = Classification), size = 2) +
    theme_bw())

(richness.ecol.ARG.rare <- plot_richness(ecol.ARG.rare, x = "SoilTemperature_to10cm", measures = "Observed") +
    geom_point(aes(shape = Classification), size = 2) +
    theme_bw())
multiplot(richness.ecol.ARG.rare, richness.ecol.rplB.rare, cols = 2)

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

#mann whitney u test for significance
#perform test of gene abundance compared to 
#soil history (Classification)
richness.mwu.classification <- subset(richness.data, Classification !="Reference") %>% group_by(GeneGroup) %>% do(tidy(wilcox.test(abs(.$Observed)~.$Classification, paired = FALSE)))

#save table
write.table(richness.mwu.classification, paste(wd, "/output/richness.1.mannwhitneyu.csv", sep = ""), row.names = FALSE, sep = ",", quote = FALSE)

#relativize rarefied datasets before beta-diversity
ecol.rplB.rareREL <-  transform_sample_counts(ecol.rplB.rare, function(x) x/sum(x))
ecol.ARG.rareREL <-  transform_sample_counts(ecol.ARG.rare, function(x) x/sum(x))

#########################################################
#B-DIVERSITY ANALYSIS OF ARGs AND AsRGs clustered at 0.01#
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

#make space into a distance matrix
space.d=as.dist(space, diag = TRUE, upper = TRUE)

#mantel tests v. space
mantel(ecol.rplB.rareREL.d,space.d, method = "spear")
mantel(ecol.ARG.rareREL.d,space.d, method = "spear")

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
gene_abundance_blast_final <- gene_abundance_blast_final[order(rowSums(dcast),decreasing=TRUE),]

#plot based on taxonomy (separate rplB)
ggplot(subset(gene_abundance_blast_final), aes(x = phylum, y = meanAbund, fill = Group, group = Group)) +
  geom_point(aes(color = Group)) +
  facet_grid(~Classification) +
  theme(axis.text.x = element_text(angle = 90, size = 8, 
                                   hjust=0.95,vjust=0.2))

gene_abundance_cast <- dcast(gene_abundance_blast_final, Classification+phylum ~ Group, fun.aggregate = sum)
la <- subset(gene_abundance_cast, phylum %in% top.phy)

gene_abundance_cast %>% group_by(phylum) %>% do(tidy(cor.test(.$rplB, .$AntibioticResistance)))

ggplot(subset(gene_abundance_cast, phylum %in% top.phy), aes(x = rplB, y = ArsenicResistance)) +
  geom_point(aes(color = phylum)) +
  geom_smooth() +
  facet_wrap(~Gene)


#### septtt
#join with BLAST data
gene_abundance_blast$Site <- gsub("Cen", "cen", gene_abundance_blast$Site)
gene_abundance_blast_wkr <- gene_abundance_blast %>%
  rename(OTU = number) %>%
  left_join(blast.final, by = "OTU") %>%
  group_by(Classification, Site, SoilTemperature_to10cm, Group, Gene.x,OTU, phylum, class) %>%
  summarise(Tot= sum(RelativeAbundance))

gene_abundance_blast_wkr2 <- gene_abundance_blast %>%
  rename(OTU = number) %>%
  left_join(blast.final, by = "OTU") %>%
  group_by(Classification, Site, SoilTemperature_to10cm, Group, Gene.x, phylum, class) %>%
  summarise(Tot= sum(RelativeAbundance))
a <- c("dfra12", "rplB")
ggplot(subset(gene_abundance_blast_wkr, Gene.x %in% a), aes(x = SoilTemperature_to10cm, y = Tot)) +
  geom_point(aes(color = Gene.x)) +
  geom_line(aes(color = Gene.x)) +
  facet_wrap(~class)

gene_abundance_blast_wkr3 <- gene_abundance_blast_wkr
gene_abundance_blast_wkr3$Site <- gsub("cen", "Cen", gene_abundance_blast_wkr3$Site)
gene_abundance_blast_wkr3$phylum <- with(gene_abundance_blast_wkr3, ifelse( phylum == "Proteobacteria", class, phylum)) 

gene_abundance_blast_wkr3 <- gene_abundance_blast_wkr3 %>%
  rename(Taxon = phylum) %>%
  left_join(rplB, by = c("Site", "Taxon"))
gene_abundance_blast_wkr3$Fraction.Abundance[is.na(gene_abundance_blast_wkr3$Fraction.Abundance)] <- 0


haha <- subset(gene_abundance_blast_wkr3, Gene.x !="rplB") %>% group_by(Gene.x, Taxon) %>% do(tidy(cor.test(.$Tot,.$Fraction.Abundance, paired = TRUE, method = "spearman")))

gene_abundance_blast_wkr3$Taxon <- ifelse(gene_abundance_blast_wkr3$class, gene_abundance_blast_wkr3$Taxon, gene_abundance_blast_wkr3$class)


                                         
gene_abundance_blast_wkr$Site <- factor(gene_abundance_blast_wkr$Site, 
                                        gene_abundance_blast_wkr$Site[order(gene_abundance_blast_wkr$SoilTemperature_to10cm)])


no <- c("rplB", "tetW", "tetX", "tetA", "AAC6-Ia", "CEP")
ggplot(subset(gene_abundance_blast_wkr, !Gene.x %in% no), aes(x = Site, y = RelativeAbundance, fill = class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Gene.x, scales = "free_y")


#correlate OTUs with Temp
spearman_otuTemp <- subset(gene_abundance_blast_wkr, Gene.x != "rplB") %>% group_by(OTU) %>% do(tidy(cor.test(.$RelativeAbundance, log(.$SoilTemperature_to10cm), method = "spearman")))

#save as table for supplemental material
write.table(spearman_richness, file = paste(wd, "/output/richness.Temp.correlations.csv", sep = ""), sep = ",", row.names = TRUE, quote = FALSE)



ggplot(subset(gene_abundance_blast_wkr), aes(x = SoilTemperature_to10cm, y = Tot, color = class)) +
  geom_smooth(aes(group = OTU),se = FALSE, method = "loess") +
  facet_wrap(~Gene.x, scales = "free_y")

ggplot(subset(gene_abundance_blast_wkr, Gene.x == "sul2"), aes(SoilTemperature_to10cm, Tot, color = phylum)) + 
  geom_point() + 
  lapply(split(gene_abundance_blast_wkr, gene_abundance_blast_wkr$phylum), function(.data){stat_function(aes(color = phylum), .data, fun = function(x){signal::pchip(.data$SoilTemperature_to10cm, .data$Tot, x)})})
