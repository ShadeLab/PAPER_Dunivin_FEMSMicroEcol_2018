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
filenames.1 <- list.files(pattern="*_rformat_dist_0.1.txt")

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
mixed.03 <- c("arxA_12", "arxA_04", "arrA_1", "arrA_2", "arrA_3", "arrA_5", "arrA_6", "arrA_8", "arrA_9", "aioA_18", "aioA_30", "aioA_36", "aioA_44", "aioA_52", "aioA_53", "aioA_57")

#remove OTUs that are gene matches
otu_table <- otu_table[,!names(otu_table) %in% mixed.03]

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

#list OTUs present in less than 2 samples
abund <- otu_table_normPA[which(rowSums(otu_table_normPA) > 3),]

#remove OTUs with presence in less than 4 sites
otu_table_norm.slim <- otu_table_norm_annotated.t[which(rownames(otu_table_norm_annotated.t) %in% rownames(abund)),]

#transpose dataset
otu_table_norm.slim.t <- t(otu_table_norm.slim)

#perform network analysis without rplB!
#remove column based on pattern (rplB)
otu_table_norm.slim.t.genes <- otu_table_norm.slim.t[, -grep("rplB", colnames(otu_table_norm.slim.t))]

#find correlations between contigs!
corr.genes <- corr.test(otu_table_norm.slim.t, 
                        method = "spearman", adjust = "fdr", alpha = 0.1)

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

r$Group[r$Gene == "rplB"] <- "Organism"
r$Shape[r$Group == "Organism"] <- "square"
r$Group[is.na(r$Group)] <- "Metadata"
r$Shape[r$Group == "Metadata"] <- "diamond"
r$Shape[r$Group == "AntibioticResistance"] <- "triangle"
shapes <- as.vector(r$Shape)

#make rplB and metadata gene colors
r$gene.color[r$Gene == "rplB"] <- "#fffac8"
r$gene.color[r$Group == "Metadata"] <- "#D3D3D3"

#save list of remaining OTUs (for BLAST)
write_lines(r$gene, "network_contigs_0.1.txt")

#read in BLAST output results
setwd(paste(wd, "/../networks/data", sep = ""))
blast.names <- list.files(pattern="*0.1.txt")
blast <- do.call(rbind, lapply(blast.names, function(X) {
  data.frame(id = basename(X), read_delim(X, delim = "\t", col_names = FALSE))}))


#tidy blast results
blast$id <- gsub("arsC_", "arsC", blast$id)
blast.tidy <- blast %>%
  separate(id, into = c("results", "Gene", "Clust"), sep = "_") %>%
  separate(X3, into = c("accno", "other"), sep = " coded_by=") %>%
  separate(other, into = c("list", "organism"), sep = ",organism=") %>%
  separate(organism, into = c("organism", "definition"), sep = ",definition=") %>%
  rename(OTU = X1, evalue = X4, percid = X5) %>%
  select(-c(results, X2))

#subset rplB results (pre-classified)
blast.rplB <- subset(blast.tidy, Gene == "rplB")
blast.all <- subset(blast.tidy, Gene != "rplB")

#get taxanomic information for each contig
#will take ~2 minutes!
blast.all[,10:14] <- tax_name(blast.all$organism, get = c("genus", "class", "phylum"), db = "ncbi")
blast.all <- blast.all %>%
  select(Gene, OTU, phylum)

#tidy rplB information (ie get phyla)
blast.rplB.tidy <- blast.rplB %>%
  separate(accno, into = c("kingdom", "phylum"), sep = ";") %>%
  select(Gene, OTU, phylum)

#join together rplB and blast data 
blast.final <- rbind(blast.rplB.tidy, blast.all)

#replace "OTU" with gene name so it can be joined with shape
blast.final$OTU <- with(blast.final, str_replace_all(blast.final$OTU, fixed("OTU"), blast.final$Gene))

#add leading zeros to r otus
r.final <- r %>%
  separate(OTU, into = c("beginning", "number"), sep = "_")

r.final$number <- sprintf("%04s", r.final$number)
r.final$OTU <- paste(r.final$beginning, r.final$number, sep = "_")

#join aesthetic (shape, taxonomy) information
aesthetics <- r.final %>%
  left_join(blast.final, by = "OTU")
aesthetics$phylum[is.na(aesthetics$phylum)] <- "unknown"
aesthetics$phylum[aesthetics$phylum == "artificialsequences"] <- "unknown"
aesthetics$phylum[aesthetics$phylum == "metagenomes"] <- "unknown"

#read in colors for phyla
phylum.colors <- read_delim(paste(wd, "/../networks/data/phylum_colors.txt", sep = ""), delim = " ", col_names = c("phy.color", "phylum"))

#add phylum colors to aesthetics
aesthetics <- aesthetics %>%
  left_join(phylum.colors, by = "phylum")


#make network of correlations
clust.network <- qgraph(corr.genes$r, minimum = "sig", sampleSize=13, 
                        layout = "spring", details = TRUE,
                        graph = "cor",label.cex = 0.5,
                        alpha = 0.01, graph = "fdr", labels = aesthetics$number,  label.scale.equal = TRUE, label.scale = FALSE,shape = aesthetics$Shape, node.resolution = 500, color = aesthetics$phy.color, negDashed = TRUE, curve = 0.2, curveAll = TRUE, overlay = FALSE)


clust.network <- qgraph(corr.genes$r, minimum = "sig", sampleSize=13, 
                         details = FALSE, layout = "spring",
                        graph = "cor",label.cex = 0.5,
                        alpha = 0.1, graph = "fdr", labels = aesthetics$number,  label.scale.equal = TRUE, label.scale = FALSE,shape = aesthetics$Shape, node.resolution = 500, negDashed = TRUE, curve = 0.2, curveAll = TRUE, overlay = TRUE, groups = groups, palette = "ggplot2")


####################################
#AsRG-ARG CLUSTER NETWORK WITH RPLB#
####################################
##make network with few genes
#subset data to only contain all of the AsRG ARG clusters
#list otus that are gene matches
matches <- c("tolC_07", "acr3_119", "arsM_0234", "rplB_0263","rplB_0736", "SoilTemperature_to10cm",
             "Ca_ppm", "Mg_ppm")

#remove non matches from OTU table
otu_table_norm.slim.t_matches <- otu_table_norm.slim.t[,colnames(otu_table_norm.slim.t) %in% matches]

#find correlations between OTUs
corr.clust <- corr.test(otu_table_norm.slim.t, 
                        method = "spearman", adjust = "fdr", alpha = 0.1)

corr.clust.r <- as.matrix(print(corr.clust$r, long = TRUE))
corr.clust.p <- as.matrix(print(corr.clust$p, long = TRUE))

corr.clust.r[which(corr.clust.p > 0.1)] <- 0

#save correlation data
write.table(corr.clust$r, paste(wd, "/output/corr_table.rplB.txt", sep = ""), quote = FALSE)

#make network of correlations
clust.network <- qgraph(corr.clust$r, minimum = "sig", sampleSize=13, 
                        layout = "spring", details = TRUE,
                        graph = "cor", label.cex = 1,
                        alpha = 0.1)

####################################################
#EXAMINE ABUNDANCE OF CO-OCCURRING RESISTANCE GENES#
####################################################

#plot interesting trends
match_trends <- otu_table_norm.slim.t_matches %>%
  melt(value.name = "Abundance") %>%
  rename(Site = Var1, OTU = Var2) %>%
  left_join(meta, by = "Site")

#order based on group
match_trends$Site <- factor(match_trends$Site, 
                            match_trends$Site[order(match_trends$SoilTemperature_to10cm)])

color1 <- c("#8DD3C7", "#FFFFB3", "#F7BD84", "#9CD1FF","#E7298A", "#FB8072")
match1 <- c("tolC_07", "acr3_119", "arsM_0234","rplB_0263","rplB_0736")
match_trends1 <- match_trends[which(match_trends$OTU %in% match1),]
(mt1 <- ggplot(match_trends1, aes(x = Site, y = Abundance, fill = OTU)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = color1) +
    theme_classic(base_size = 12) +
    ylim(0,0.06) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=0.95,vjust=0.9)))

color2 <- c("#8DD3C7", "#FFFFB3","#9CD1FF", "#7570B3")
match2 <- c("rplB_0692", "acr3_002", "arsM_296", "dfra12_038")
match_trends2 <- match_trends[which(match_trends$OTU %in% match2),]
(mt2 <- ggplot(match_trends2, 
               aes(x = Site, y = Abundance, fill = OTU)) +
    geom_bar(stat = "identity", color = "black")  +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = color2) +
    ylim(0,0.06) +
    theme(axis.text.x = element_text(angle = 45,  
                                     hjust=0.95,vjust=0.9)))

#save abundance of ARG-AsRG clusters
ggsave(mt1, filename = paste(wd, "/figures/cluster.abundance.grp1.eps", sep = ""), width = 5, height = 3, units = "in")
ggsave(mt2, filename = paste(wd, "/figures/cluster.abundance.grp2.eps", sep = ""), width = 5, height = 3, units = "in")

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

########################################################
#ECOLOGICAL ANALYSIS OF ARGs AND AsRGs clustered at 0.1#
########################################################


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
class.033 <- c("#FFFF00", "darkgreen", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000")

#check sampling depth of each matrix
rarecurve(ecol.rplB, step=1, label = FALSE, col = class.033)
rarecurve(ecol.AsRG, step=1, label = FALSE, col = class.033)
rarecurve(ecol.ARG, step=1, label = FALSE, col = class.033)

#convert contig tables into phyloseq objects
ecol.rplB.phyloseq <- otu_table(ecol.rplB, taxa_are_rows = FALSE)
ecol.AsRG.phyloseq <- otu_table(ecol.AsRG, taxa_are_rows = FALSE)
ecol.ARG.phyloseq <- otu_table(ecol.ARG, taxa_are_rows = FALSE)

#rarefy to even sampling depth 
ecol.rplB.rare <- rarefy_even_depth(ecol.rplB.phyloseq, rngseed = TRUE)
rarecurve(ecol.rplB.rare, step=1, label = FALSE, col = class.033)

ecol.ARG.rare <- rarefy_even_depth(ecol.ARG.phyloseq, rngseed = TRUE)
rarecurve(ecol.ARG.rare, step=1, label = FALSE, col = class.033)

ecol.AsRG.rare <- rarefy_even_depth(ecol.AsRG.phyloseq, rngseed = TRUE)
rarecurve(ecol.AsRG.rare, step=1, label = FALSE, col = class.033)

#calculate evenness
plieou.rplB <- data.frame(group = "rplB", Site = rownames(ecol.rplB.rare), plieou = vegan::diversity(ecol.rplB.rare, index = "shannon")/log(specnumber(ecol.rplB.rare)))

plieou.ARG <- data.frame(group = "ARG", Site = rownames(ecol.ARG.rare), plieou = vegan::diversity(ecol.ARG.rare, index = "shannon")/log(specnumber(ecol.ARG.rare)))

plieou.AsRG <- data.frame(group = "AsRG", Site = rownames(ecol.AsRG.rare), plieou = vegan::diversity(ecol.AsRG.rare, index = "shannon")/log(specnumber(ecol.AsRG.rare)))

#join all evenness information and add metadata
plieou.full <- rbind(plieou.ARG, plieou.AsRG, plieou.rplB)
plieou.full <- left_join(plieou.full, meta, by = "Site")

#plot evenness
(plieou.plot <- ggplot(plieou.full, aes(x = Classification, y = plieou)) +
    geom_boxplot() +
    ylab(label = "Evenness") +
    theme_bw() +
    facet_wrap(~group))

#test differences in evenness
#mann whitney u test for significance
#perform test of gene abundance compared to 
#soil history (Classification)
evenness.mwu.classification <- subset(plieou.full, Classification !="Reference") %>% group_by(group) %>% do(tidy(wilcox.test(abs(.$plieou)~.$Classification, paired = FALSE)))

#save table
write.table(evenness.mwu.classification, paste(wd, "/output/evenness.03.mannwhitneyu.csv", sep = ""), row.names = FALSE, sep = ",", quote = FALSE)

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

richness0.1 <- multiplot(richness.ecol.rplB.rare, richness.ecol.ARG.rare, richness.ecol.AsRG.rare, cols=3)

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
write.table(richness.mwu.classification, paste(wd, "/output/richness.03.mannwhitneyu.csv", sep = ""), row.names = FALSE, sep = ",", quote = FALSE)

#relativize rarefied datasets before beta-diversity
ecol.rplB.rareREL <-  transform_sample_counts(ecol.rplB.rare, function(x) x/sum(x))
ecol.ARG.rareREL <-  transform_sample_counts(ecol.ARG.rare, function(x) x/sum(x))
ecol.AsRG.rareREL <- transform_sample_counts(ecol.AsRG.rare, function(x) x/sum(x))

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

#########################################################
#ECOLOGICAL ANALYSIS OF ARGs AND AsRGs clustered at 0.1#
#########################################################

#make an OTU table with 0.1 cluster sizes

#temporarily change working directories
setwd(paste(wd, "/data", sep = ""))

#list filenames of interest
filenames.03 <- list.files(pattern="*_rformat_dist_0.1.txt")

#move back up directories
setwd("../..")

#make dataframes of all OTU tables
for(i in filenames.03){
  filepath <- file.path(paste(wd, "/data", sep = ""),paste(i,sep=""))
  assign(gsub("_rformat_dist_0.1.txt", ".03", i), read.delim(filepath,sep = "\t"))
}

#change OTU to gene name
colnames(acr3.03) <- naming(acr3.03)
colnames(aioA.03) <- naming(aioA.03)
colnames(arsB.03) <- naming(arsB.03)
colnames(`AAC6-Ia.03`) <- naming(`AAC6-Ia.03`)
colnames(adeB.03) <- naming(adeB.03)
colnames(arrA.03) <- naming(arrA.03)
colnames(arsA.03) <- naming(arsA.03)
colnames(arsC_glut.03) <- naming(arsC_glut.03)
colnames(arsC_thio.03) <- naming(arsC_thio.03)
colnames(arsD.03) <- naming(arsD.03)
colnames(arsM.03) <- naming(arsM.03)
colnames(arxA.03) <- naming(arxA.03)
colnames(CEP.03) <- naming(CEP.03)
colnames(ClassA.03) <- naming(ClassA.03)
colnames(ClassB.03) <- naming(ClassB.03)
colnames(ClassC.03) <- naming(ClassC.03)
colnames(dfra12.03) <- naming(dfra12.03)
colnames(rplB.03) <- naming(rplB.03)
colnames(intI.03) <- naming(intI.03)
colnames(sul2.03) <- naming(sul2.03)
colnames(tetA.03) <- naming(tetA.03)
colnames(tetW.03) <- naming(tetW.03)
colnames(tetX.03) <- naming(tetX.03)
colnames(tolC.03) <- naming(tolC.03)
colnames(vanA.03) <- naming(vanA.03)
colnames(vanH.03) <- naming(vanH.03)
colnames(vanX.03) <- naming(vanX.03)
colnames(vanZ.03) <- naming(vanZ.03)

#join together all files
otu_table.03 <- acr3.03 %>%
  left_join(aioA.03, by = "X") %>%
  left_join(`AAC6-Ia.03`, by = "X") %>%
  left_join(arrA.03, by = "X") %>%
  left_join(arsA.03, by = "X") %>%
  left_join(arsB.03, by = "X") %>%
  left_join(arsC_glut.03, by = "X") %>%
  left_join(arsC_thio.03, by = "X") %>%
  left_join(arsD.03, by = "X") %>%
  left_join(arsM.03, by = "X") %>%
  left_join(arxA.03, by = "X") %>%
  left_join(CEP.03, by = "X") %>%
  left_join(ClassA.03, by = "X") %>%
  left_join(ClassB.03, by = "X") %>%
  left_join(ClassC.03, by = "X") %>%
  left_join(dfra12.03, by = "X") %>%
  left_join(intI.03, by = "X") %>%
  left_join(rplB.03, by = "X") %>%
  left_join(sul2.03, by = "X") %>%
  left_join(tetA.03, by = "X") %>%
  left_join(tetW.03, by = "X") %>%
  left_join(tetX.03, by = "X") %>%
  left_join(tolC.03, by = "X") %>%
  left_join(vanA.03, by = "X") %>%
  left_join(vanH.03, by = "X") %>%
  left_join(vanX.03, by = "X") %>%
  left_join(vanZ.03, by = "X") %>%
  left_join(adeB.03, by = "X") %>%
  rename(Site =X) 

#list otus that are gene matches
#mixed.03 <- c("arxA_12", "arxA_04", "arrA_1", "arrA_2", "arrA_3", "arrA_5", "arrA_6", "arrA_8", "arrA_9", "aioA_18", "aioA_30", "aioA_36", "aioA_44", "aioA_52", "aioA_53", "aioA_57")

#remove OTUs that are gene matches
#otu_table <- otu_table[,!names(otu_table) %in% mixed.03]

#rename arxA column that is actually arrA (with no match)
#otu_table <- otu_table %>%
#  rename(arrA_10 = arxA_10)

#replace all NAs (from join) with zeros
otu_table.03[is.na(otu_table.03)] <- 0

#separate original (non-normalized) contig table into
#specific gene groups (rplB, ARG, AsRG)
ecol.rplB <- otu_table.03[, grep("rplB", colnames(otu_table.03))]
ecol.AsRG <- otu_table.03[, grep("ars|aio|arx|arr|acr", 
                                 colnames(otu_table.03))]
ecol.ARG <- otu_table.03[,grep("tolC|dfra|Class|tet|van|CEP|AAC|ade|sul", 
                               colnames(otu_table.03))]

#add site names as row names to each contig table 
rownames(ecol.rplB) <- otu_table[,1]
rownames(ecol.ARG) <- otu_table[,1]
rownames(ecol.AsRG) <- otu_table[,1]

#make list of 13 colors (based on classification)
class.033 <- c("#FFFF00", "darkgreen", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000")

#check sampling depth of each matrix
rarecurve(ecol.rplB, step=1, label = FALSE, col = class.033)
rarecurve(ecol.AsRG, step=1, label = FALSE, col = class.033)
rarecurve(ecol.ARG, step=1, label = FALSE, col = class.033)

#convert contig tables into phyloseq objects
ecol.rplB.phyloseq <- otu_table(ecol.rplB, taxa_are_rows = FALSE)
ecol.AsRG.phyloseq <- otu_table(ecol.AsRG, taxa_are_rows = FALSE)
ecol.ARG.phyloseq <- otu_table(ecol.ARG, taxa_are_rows = FALSE)

#rarefy to even sampling depth 
ecol.rplB.rare <- rarefy_even_depth(ecol.rplB.phyloseq, rngseed = TRUE)
rarecurve(ecol.rplB.rare, step=1, label = FALSE, col = class.033)

ecol.ARG.rare <- rarefy_even_depth(ecol.ARG.phyloseq, rngseed = TRUE)
rarecurve(ecol.ARG.rare, step=1, label = FALSE, col = class.033)

ecol.AsRG.rare <- rarefy_even_depth(ecol.AsRG.phyloseq, rngseed = TRUE)
rarecurve(ecol.AsRG.rare, step=1, label = FALSE, col = class.033)

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
    ylim(30,170) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

(richness.ecol.ARG.rare <- plot_richness(ecol.ARG.rare, x = "Classification", measures = "Observed") +
    geom_boxplot() +
    theme_bw() +
    ylim(30,170) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

(richness.ecol.AsRG.rare <- plot_richness(ecol.AsRG.rare, x = "Classification", measures = "Observed") +
    geom_boxplot() +
    theme_bw() +
    ylim(30,170) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

richness0.1 <- multiplot(richness.ecol.rplB.rare, richness.ecol.ARG.rare, richness.ecol.AsRG.rare, cols=3)

ggsave(richness0.1, filename = paste(wd, "/figures/richness.png",sep = ""), width = 6, units = "in")

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
richness.mwu.classification <- subset(richness.data, Classification !="Reference") %>% group_by(GeneGroup) %>% do(tidy(t.test(abs(.$Observed)~.$Classification, paired = FALSE)))

#save table
write.table(richness.mwu.classification, paste(wd, "/output/richness.03.mannwhitneyu.csv", sep = ""), row.names = FALSE, sep = ",", quote = FALSE)

#relativize rarefied datasets before beta-diversity
ecol.rplB.rareREL <-  transform_sample_counts(ecol.rplB.rare, function(x) x/sum(x))
ecol.ARG.rareREL <-  transform_sample_counts(ecol.ARG.rare, function(x) x/sum(x))
ecol.AsRG.rareREL <- transform_sample_counts(ecol.AsRG.rare, function(x) x/sum(x))

#plot Bray Curtis ordination for rplB
ord.rplB.bray <- ordinate(ecol.rplB.rareREL, method="PCoA", distance="bray")
(bc.ord.rplB=plot_ordination(ecol.rplB.rareREL, ord.rplB.bray, color="SoilTemperature_to10cm",
                             title="Bray Curtis", shape = "Classification") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    geom_point(size=5, alpha = 0.5) +
    theme_light(base_size = 12))

#plot Bray Curtis ordination for ARGs
ord.ARG.bray <- ordinate(ecol.ARG.rareREL, method="PCoA", distance="bray")
(bc.ord.ARG=plot_ordination(ecol.ARG.rareREL, ord.ARG.bray, color="SoilTemperature_to10cm",
                            title="Bray Curtis", shape = "Classification") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    geom_point(size=5, alpha = 0.5) +
    theme_light(base_size = 12))

#plot Bray Curtis ordination for AsRGs
ord.AsRG.bray <- ordinate(ecol.AsRG.rareREL, method="PCoA", distance="bray")
(bc.ord.AsRG=plot_ordination(ecol.AsRG.rareREL, ord.AsRG.bray, color="SoilTemperature_to10cm",
                             title="Bray Curtis", shape = "Classification") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    geom_point(size=5, alpha = 0.5) +
    theme_light(base_size = 12))

bray.curtis <- multiplot(bc.ord.rplB, bc.ord.ARG, bc.ord.AsRG, ncol = 1)

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

#mantel tests v. space
mantel(ecol.rplB.rareREL.d,space.d, method = "spear")
mantel(ecol.ARG.rareREL.d,space.d, method = "spear")
mantel(ecol.AsRG.rareREL.d,space.d, method = "spear")

#####################################
#0.1 CLUSTERING ECOLOGICAL ANALYSIS#
#####################################

#make an OTU table with 0.1 cluster sizes

#temporarily change working directories
setwd(paste(wd, "/data", sep = ""))

#list filenames of interest
filenames.03 <- list.files(pattern="*_rformat_dist_0.1.txt")

#move back up directories
setwd("../..")

#make dataframes of all OTU tables
for(i in filenames.03){
  filepath <- file.path(paste(wd, "/data", sep = ""),paste(i,sep=""))
  assign(gsub("_rformat_dist_0.1.txt", ".03", i), read.delim(filepath,sep = "\t"))
}

#change OTU to gene name
colnames(acr3.03) <- naming(acr3.03)
colnames(aioA.03) <- naming(aioA.03)
colnames(arsB.03) <- naming(arsB.03)
colnames(`AAC6-Ia.03`) <- naming(`AAC6-Ia.03`)
colnames(adeB.03) <- naming(adeB.03)
colnames(arrA.03) <- naming(arrA.03)
colnames(arsA.03) <- naming(arsA.03)
colnames(arsC_glut.03) <- naming(arsC_glut.03)
colnames(arsC_thio.03) <- naming(arsC_thio.03)
colnames(arsD.03) <- naming(arsD.03)
colnames(arsM.03) <- naming(arsM.03)
colnames(arxA.03) <- naming(arxA.03)
colnames(CEP.03) <- naming(CEP.03)
colnames(ClassA.03) <- naming(ClassA.03)
colnames(ClassB.03) <- naming(ClassB.03)
colnames(ClassC.03) <- naming(ClassC.03)
colnames(dfra12.03) <- naming(dfra12.03)
colnames(rplB.03) <- naming(rplB.03)
colnames(intI.03) <- naming(intI.03)
colnames(sul2.03) <- naming(sul2.03)
colnames(tetA.03) <- naming(tetA.03)
colnames(tetW.03) <- naming(tetW.03)
colnames(tetX.03) <- naming(tetX.03)
colnames(tolC.03) <- naming(tolC.03)
colnames(vanA.03) <- naming(vanA.03)
colnames(vanH.03) <- naming(vanH.03)
colnames(vanX.03) <- naming(vanX.03)
colnames(vanZ.03) <- naming(vanZ.03)

#join together all files
otu_table.03 <- acr3.03 %>%
  left_join(aioA.03, by = "X") %>%
  left_join(`AAC6-Ia.03`, by = "X") %>%
  left_join(arrA.03, by = "X") %>%
  left_join(arsA.03, by = "X") %>%
  left_join(arsB.03, by = "X") %>%
  left_join(arsC_glut.03, by = "X") %>%
  left_join(arsC_thio.03, by = "X") %>%
  left_join(arsD.03, by = "X") %>%
  left_join(arsM.03, by = "X") %>%
  left_join(arxA.03, by = "X") %>%
  left_join(CEP.03, by = "X") %>%
  left_join(ClassA.03, by = "X") %>%
  left_join(ClassB.03, by = "X") %>%
  left_join(ClassC.03, by = "X") %>%
  left_join(dfra12.03, by = "X") %>%
  left_join(intI.03, by = "X") %>%
  left_join(rplB.03, by = "X") %>%
  left_join(sul2.03, by = "X") %>%
  left_join(tetA.03, by = "X") %>%
  left_join(tetW.03, by = "X") %>%
  left_join(tetX.03, by = "X") %>%
  left_join(tolC.03, by = "X") %>%
  left_join(vanA.03, by = "X") %>%
  left_join(vanH.03, by = "X") %>%
  left_join(vanX.03, by = "X") %>%
  left_join(vanZ.03, by = "X") %>%
  left_join(adeB.03, by = "X") %>%
  rename(Site =X) 

#list otus that are gene matches
#mixed.03 <- c("arxA_12", "arxA_04", "arrA_1", "arrA_2", "arrA_3", "arrA_5", "arrA_6", "arrA_8", "arrA_9", "aioA_18", "aioA_30", "aioA_36", "aioA_44", "aioA_52", "aioA_53", "aioA_57")

#remove OTUs that are gene matches
#otu_table <- otu_table[,!names(otu_table) %in% mixed.03]

#rename arxA column that is actually arrA (with no match)
#otu_table <- otu_table %>%
#  rename(arrA_10 = arxA_10)

#replace all NAs (from join) with zeros
otu_table.03[is.na(otu_table.03)] <- 0

#separate original (non-normalized) contig table into
#specific gene groups (rplB, ARG, AsRG)
ecol.rplB <- otu_table.03[, grep("rplB", colnames(otu_table.03))]
ecol.AsRG <- otu_table.03[, grep("ars|aio|arx|arr|acr", 
                                 colnames(otu_table.03))]
ecol.ARG <- otu_table.03[,grep("tolC|dfra|Class|tet|van|CEP|AAC|ade|sul", 
                               colnames(otu_table.03))]

#add site names as row names to each contig table 
rownames(ecol.rplB) <- otu_table[,1]
rownames(ecol.ARG) <- otu_table[,1]
rownames(ecol.AsRG) <- otu_table[,1]

#make list of 13 colors (based on classification)
class.033 <- c("#FFFF00", "darkgreen", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000", "#FF0000", "#FFFF00", "#FF0000")

#check sampling depth of each matrix
rarecurve(ecol.rplB, step=1, label = TRUE, col = class.033)
rarecurve(ecol.AsRG, step=1, label = TRUE, col = class.033)
rarecurve(ecol.ARG, step=1, label = FALSE, col = class.033)

#convert contig tables into phyloseq objects
ecol.rplB.phyloseq <- otu_table(ecol.rplB, taxa_are_rows = FALSE)
ecol.AsRG.phyloseq <- otu_table(ecol.AsRG, taxa_are_rows = FALSE)
ecol.ARG.phyloseq <- otu_table(ecol.ARG, taxa_are_rows = FALSE)

#rarefy to even sampling depth 
ecol.rplB.rare <- rarefy_even_depth(ecol.rplB.phyloseq, rngseed = TRUE)
rarecurve(ecol.rplB.rare, step=1, label = FALSE, col = class.033)

ecol.ARG.rare <- rarefy_even_depth(ecol.ARG.phyloseq, rngseed = TRUE)
rarecurve(ecol.ARG.rare, step=1, label = FALSE, col = class.033)

ecol.AsRG.rare <- rarefy_even_depth(ecol.AsRG.phyloseq, rngseed = TRUE)
rarecurve(ecol.AsRG.rare, step=1, label = FALSE, col = class.033)

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
    ylim(30, 230) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

(richness.ecol.ARG.rare <- plot_richness(ecol.ARG.rare, x = "Classification", measures = "Observed") +
    geom_boxplot() +
    theme_bw() +
    ylim(30, 230) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

(richness.ecol.AsRG.rare <- plot_richness(ecol.AsRG.rare, x = "Classification", measures = "Observed") +
    geom_boxplot() +
    theme_bw() +
    ylim(30, 230) +
    theme(axis.text.x = element_text(angle = 45, size = 14, 
                                     hjust=0.95)))

richness <- multiplot(richness.ecol.rplB.rare, richness.ecol.ARG.rare, richness.ecol.AsRG.rare, cols=3)

ggsave(richness, filename = paste(wd, "/figures/richness.png",sep = ""), width = 6, units = "in")

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
write.table(richness.mwu.classification, paste(wd, "/output/richness.03.mannwhitneyu.csv", sep = ""), row.names = FALSE, sep = ",", quote = FALSE)


#relativize rarefied datasets before beta-diversity
ecol.rplB.rareREL <-  transform_sample_counts(ecol.rplB.rare, function(x) x/sum(x))
ecol.ARG.rareREL <-  transform_sample_counts(ecol.ARG.rare, function(x) x/sum(x))
ecol.AsRG.rareREL <- transform_sample_counts(ecol.AsRG.rare, function(x) x/sum(x))

#plot Bray Curtis ordination for rplB
ord.rplB.bray <- ordinate(ecol.rplB.rareREL, method="PCoA", distance="bray")
(bc.ord.rplB=plot_ordination(ecol.rplB.rareREL, ord.rplB.bray, color="SoilTemperature_to10cm",
                             title="Bray Curtis", shape = "Classification") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                          guide_legend(title="Temperature (°C)")) +
    geom_point(size=5, alpha = 0.5) +
    theme_light(base_size = 12))

#plot Bray Curtis ordination for ARGs
ord.ARG.bray <- ordinate(ecol.ARG.rareREL, method="PCoA", distance="bray")
(bc.ord.ARG=plot_ordination(ecol.ARG.rareREL, ord.ARG.bray, color="SoilTemperature_to10cm",
                            title="Bray Curtis", shape = "Classification") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide = FALSE) +
    geom_point(size=5, alpha = 0.5) +
    theme_light(base_size = 12))

#plot Bray Curtis ordination for AsRGs
ord.AsRG.bray <- ordinate(ecol.AsRG.rareREL, method="PCoA", distance="bray")
(bc.ord.AsRG=plot_ordination(ecol.AsRG.rareREL, ord.AsRG.bray, color="SoilTemperature_to10cm",
                             title="Bray Curtis", shape = "Classification") +
    scale_color_gradientn(colours=GnYlOrRd(5), guide = FALSE) +
    geom_point(size=5, alpha = 0.5) +
    theme_light(base_size = 12))

bray.curtis <- multiplot(bc.ord.rplB, bc.ord.ARG, bc.ord.AsRG, ncol = 1)

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

#mantel tests v. space
mantel(ecol.rplB.rareREL.d,space.d, method = "spear")
mantel(ecol.ARG.rareREL.d,space.d, method = "spear")
mantel(ecol.AsRG.rareREL.d,space.d, method = "spear")
#we do not see different results with 0.1 clustering

###?##?#?#?#?#?#?#?#?#?#?#?#?#?





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
    theme_classic(base_size = 8) +
    ylab("Gene per rplB (%)") +
    facet_wrap(~ Gene, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=0.95,vjust=0.2)))

(gene.bar.census <- ggplot(subset(data.phylum, Gene == "arsM"),
                           aes(x = Site,  y = Phylum.count*100, fill = phylum)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    theme_classic(base_size = 8) +
    ylab("Gene per rplB (%)") +
    facet_wrap(~ Gene, scales = "free_y") +
    scale_fill_manual(values = color) +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=0.95,vjust=0.2)))

(gene.bar.census <- ggplot(subset(data.phylum, Gene == "arsM"),
                           aes(x = Site,  y = Phylum.count*100, fill = phylum)) +
    geom_bar(stat = "identity") +
    theme_classic(base_size = 8) +
    ylab("Gene per rplB (%)") +
    facet_wrap(~ Gene, scales = "free_y") +
    scale_fill_manual(values = color) +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=0.95,vjust=0.2)))

#look at class level differneces
data.family <- identifiers_ncbi %>%
  rename(Temp = SoilTemperature_to10cm) %>%
  group_by(Group, Description, Gene, family, Classification, Site, Temp) %>%
  summarise(Phylum.count = sum(RelativeAbundance))

(gene.bar.census <- ggplot(subset(data.class, Gene == "arsM"),
                           aes(x = Site,  y = Phylum.count*100, fill = class)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    theme_classic(base_size = 8) +
    ylab("Gene per rplB (%)") +
    facet_wrap(~ Gene, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=0.95,vjust=0.2)))















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
    theme_classic(base_size = 8) +
    ylab("Gene per rplB (%)") +
    facet_wrap(~ Gene, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=0.95,vjust=0.2)))

(gene.bar.census <- ggplot(subset(data.phylum, Gene == "arsM"),
                           aes(x = Site,  y = Phylum.count*100, fill = phylum)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    theme_classic(base_size = 8) +
    ylab("Gene per rplB (%)") +
    facet_wrap(~ Gene, scales = "free_y") +
    scale_fill_manual(values = color) +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=0.95,vjust=0.2)))

(gene.bar.census <- ggplot(subset(data.phylum, Gene == "arsM"),
                           aes(x = Site,  y = Phylum.count*100, fill = phylum)) +
    geom_bar(stat = "identity") +
    theme_classic(base_size = 8) +
    ylab("Gene per rplB (%)") +
    facet_wrap(~ Gene, scales = "free_y") +
    scale_fill_manual(values = color) +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=0.95,vjust=0.2)))

#look at class level differneces
data.family <- identifiers_ncbi %>%
  rename(Temp = SoilTemperature_to10cm) %>%
  group_by(Group, Description, Gene, family, Classification, Site, Temp) %>%
  summarise(Phylum.count = sum(RelativeAbundance))

(gene.bar.census <- ggplot(subset(data.class, Gene == "arsM"),
                           aes(x = Site,  y = Phylum.count*100, fill = class)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    theme_classic(base_size = 8) +
    ylab("Gene per rplB (%)") +
    facet_wrap(~ Gene, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=0.95,vjust=0.2)))



