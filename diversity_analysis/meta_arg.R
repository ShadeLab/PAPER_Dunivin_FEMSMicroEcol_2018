##########################
#READ IN DATA, SET UP ENV#
##########################

#read dependencies
library(phyloseq)
library(vegan)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(taxize)

#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis
wd <- print(getwd())

#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_JGI_map.txt", 
                                    sep=""), sep=" ", header=TRUE))

#read in microbe census data
census <- read_delim(file = paste(wd, "/data/microbe_census.txt", sep = ""),
                     delim = "\t", col_types = list(col_character(),
                                                    col_number(),
                                                    col_number(), 
                                                    col_number()))

#read in gene classification data
gene <- read_delim(paste(wd, "/data/gene_classification.txt",  sep=""), 
                   delim = "\t", col_names = TRUE)

#make color pallette for Centralia temperatures
GnYlOrRd <- colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

####################################
#READ IN AND SET UP DATA#
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

####################################
#EXAMINE rplB ACROSS CHRONOSEQUENCE#
####################################

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

#save summarised data for future analyses
write.table(x = summarised.rplB, file = paste(wd, "/output/rplB.summary.scg.txt", sep = ""), row.names = FALSE)

#decast for abundance check
dcast <- acast(rplB, Taxon ~ Site, value.var = "Fraction.Abundance")

#call na's zeros
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
                theme(axis.text.x = element_text(angle = 90, size = 12, 
                                                 hjust=0.95,vjust=0.2))))

#save plot
ggsave(phylum.plot, filename = paste(wd, "/figures/phylum.responses.png", sep=""), 
       width = 5, height = 5)

###################################################
#EXAMINE As resistance genes ACROSS CHRONOSEQUENCE#
###################################################

#Tidy gene data
data.tidy <- data %>%
  separate(col = Taxon, into = c("Code", "Organism"), sep = "organism=") %>%
  separate(col = Organism, into = c("Organism", "Definition"), sep = ",definition=") %>%
  select(Site, Gene, Organism:Fraction.Abundance) %>%
  group_by(Gene, Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.tidy$Fraction.Abundance <- as.numeric(data.tidy$Fraction.Abundance)
data.tidy$Abundance <- as.numeric(data.tidy$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.total <- data.tidy %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#make column for organism name and join with microbe census data and normalize to it
data.annotated <- data.tidy %>%
  left_join(census, by = "Site") %>%
  left_join(summarised.rplB, by = "Site") %>%
  left_join(meta, by = "Site") %>%
  left_join(gene, by = "Gene") %>%
  rename(Temp = SoilTemperature_to10cm) %>%
  select(Site:GE, rplB, As_ppm, Temp, Classification:Group) %>%
  mutate(Normalized.Abundance.census = Abundance / GE, 
         Normalized.Abundance.rplB = Abundance / rplB)

#summarise data to get number of genes per gene per site
data.site <- data.annotated %>%
  group_by(Gene, Site) %>%
  summarise(Count = sum(Abundance), 
            Count.rplB = sum(Normalized.Abundance.rplB),
            Count.census = sum(Normalized.Abundance.census))

#add taxanomic information 
#data.ncbi <- tax_name(query = data.annotated$Organism, 
#                      get = c("genus", "class", "phylum"), db = "ncbi")


#label query "Organism" for joining purposes
#data.ncbi$Organism <- data.ncbi$query

#save this table since the above step takes a long time
#write.table(data.ncbi, file = paste(wd, "/output/ncbi.taxonomy.txt", sep = ""), 
#            row.names = FALSE)

#read in NCBI data
data.ncbi <- read_delim(paste(wd, "/output/ncbi.taxonomy.txt", sep = ""), delim = " ")

#join ncbi information with annotated data
#output should have same number of rows 
data.annotated.ncbi <- data.annotated %>%
  left_join(data.ncbi, by = "Organism") %>%
  unique()

#replace NA in phylum with unknown
data.annotated.ncbi$phylum[is.na(data.annotated.ncbi$phylum)] = "Unknown"

#call NA class by phyla
data.annotated.ncbi$class[is.na(data.annotated.ncbi$class)] <- as.character(data.annotated.ncbi$phylum[is.na(data.annotated.ncbi$class)])

#call NA genus by class (may be phyla in cases where class was NA)
data.annotated.ncbi$genus[is.na(data.annotated.ncbi$genus)] <- as.character(data.annotated.ncbi$class[is.na(data.annotated.ncbi$genus)])

#order based on temperature
data.annotated.ncbi$Site <- factor(data.annotated.ncbi$Site, 
                           levels = data.annotated.ncbi$Site[order(data.annotated.ncbi$Temp)])

##############################
#EXAMINE PHYLUM LEVEL CHANGES#
##############################

#look at phylum level differneces
data.phylum <- data.annotated.ncbi %>%
  group_by(Group, Description, Gene, phylum, Classification, Site, Temp) %>%
  summarise(Phylum.count = sum(Normalized.Abundance.census))

#prep colors for phylum diversity
color <- c("#FF7F00", "#7570B3", "#CAB2D6", "#FBB4AE", "#F0027F", "#BEBADA", "#E78AC3", "#A6D854", "#B3B3B3", "#386CB0", "#BC80BD", "#FFFFCC", "#BF5B17", "#984EA3", "#CCCCCC", "#FFFF99", "#B15928", "#F781BF", "#FDC086", "#A6CEE3", "#FDB462", "#FED9A6", "#E6AB02", "#E31A1C", "#B2DF8A", "#377EB8", "#FCCDE5", "#80B1D3", "#FFD92F", "#33A02C", "#66C2A5", "#666666", "black", "brown")

#order genes by group
data.phylum$Gene <- factor(data.phylum$Gene, 
                                   levels = data.phylum$Gene[order(data.phylum$Group)])

#plot 
(gene.bar.census <- ggplot(data.phylum, aes(x = Site, 
                                            y = Phylum.count*100, fill = phylum)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = color) +
    theme_classic(base_size = 25) +
    ylab("Genome Equivalents with Gene (%)") +
    facet_wrap(~ Gene) +
    theme(axis.text.x = element_text(angle = 90, size = 15, hjust=0.95,vjust=0.2)))

#save plot
ggsave(gene.bar.census, filename = paste(wd, "/figures/phylum.abundance.by_gene.png", sep=""), height = 20, width = 25)

#examine antibiotic v arsenic resistance genes in Centralia (ie remove intI)
data.phylum.ni <- data.phylum[-which(data.phylum$Gene == "intI"),]

#plot grouped by gene functional group 
(arg.asrg.bar <- ggplot(data.phylum.ni, aes(x = Site, 
                                            y = Phylum.count*100, fill = Gene)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = color) +
    theme_classic(base_size = 15) +
    ylab("Percent of Genome Equivalents") +
    facet_wrap(~Group) +
    theme(axis.text.x = element_text(angle = 90, size = 15, hjust=0.95,vjust=0.2)))

#save plot
ggsave(arg.asrg.bar, filename = paste(wd, "/figures/arg_asrg.abundance.by_gene.png", sep=""), height = 6)

(func.bar <- ggplot(data.phylum, aes(x = Site, 
                                            y = Phylum.count*100, fill = Gene)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = color) +
    theme_classic(base_size = 12) +
    ylab("Percent of Genome Equivalents") +
    facet_wrap( ~ Description,  scales = "free_y", ncol = 2) +
    theme(axis.text.x = element_text(angle = 90, size = 12, hjust=0.95,vjust=0.2)))

#save plot
ggsave(func.bar, filename = paste(wd, "/figures/funct.abundance.by_gene.png", 
                                  sep=""), width = 9)

#summarise data by classification
data.class <- data.phylum %>%
  ungroup() %>%
  group_by(Site, Temp, Classification, Group, Description, Gene) %>%
  summarise(Count = sum(Phylum.count))

#remove rows with v. low sample genes
data.phylum.class <- data.phylum.class[-which(data.phylum.class$Gene == "tetA"),]
data.phylum.class <- data.phylum.class[-which(data.phylum.class$Gene == "tetW"),]
data.phylum.class <- data.phylum.class[-which(data.phylum.class$Gene == "tetX"),]
data.phylum.class <- data.phylum.class[-which(data.phylum.class$Gene == "vanT"),]
data.phylum.class <- data.phylum.class[-which(data.phylum.class$Gene == "arsB"),]
data.phylum.class <- data.phylum.class[-which(data.phylum.class$Gene == "CAT"),]
data.phylum.class <- data.phylum.class[-which(data.phylum.class$Gene == "AAC3-Ia"),]
data.phylum.class <- data.phylum.class[-which(data.phylum.class$Gene == "AAC6-Ia"),]
data.phylum.class <- data.phylum.class[-which(data.phylum.class$Gene == "AAC6-II"),]

#make boxplot of gene abundance in different soils 
(boxplot.genes <- ggplot(data.phylum.class, aes(x = Classification, 
                                                y = Count*100)) +
  geom_boxplot() +
  geom_jitter(aes(color = Temp), size = 2) +
  facet_wrap( ~ Gene,  scales = "free_y") +
  ylab("Genome Equivalents with Gene (%)") +
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                        guide_legend(title="Temperature (°C)")) +
  theme_classic(base_size = 12))

#save plot
ggsave(boxplot.genes, filename = paste(wd, "/figures/boxplot.by_gene.png", 
                                  sep=""), height = 6, width = 12.5)

######################
#CORRELATION ANALYSES#
######################

#order metadata based on temp to match matrix
data.phylum.class$Site <- factor(data.phylum.class$Site, levels = data.phylum.class$Site[order(meta$Site)])

#decast for abundance check but first make Site a character
data.phylum.class$Site <- as.character(data.phylum.class$Site)
cast.gene <- data.frame(acast(data.phylum.class, Site ~ Gene, value.var = "Count"))

#call na's zeros
cast.gene[is.na(cast.gene)] =0

#correlate data with temp
corr <- print(corr.test(x = cast.gene, y = meta[,c(2,16, 20:30)], method = "spearman", adjust = "fdr"), short = FALSE)
#arsenic and antibiotic resistance genes are not correlated with temperature!

######################
#MANN WHITNEY U TESTS#
######################
#remove C17 (testing recovered v reference)
data.mann <- data.phylum.class[-which(data.phylum.class$Classification == "Reference"),]

#subset data for each gene to test:
#acr3
acr3 <- subset(x = data.mann, subset = Gene == "acr3")
acr3.cast <- print(wilcox.test(acr3$Count~acr3$Classification))

#aioA
aioA <- subset(x = data.mann, subset = Gene == "aioA")
aioA.cast <- print(wilcox.test(aioA$Count~aioA$Classification))

#arsM
arsM <- subset(x = data.mann, subset = Gene == "arsM")
arsM.cast <- print(wilcox.test(arsM$Count~arsM$Classification))

#arsCglut
arsCglut <- subset(x = data.mann, subset = Gene == "arsCglut")
arsCglut.cast <- print(wilcox.test(arsCglut$Count~arsCglut$Classification))

#arsCthio
arsCthio <- subset(x = data.mann, subset = Gene == "arsCthio")
arsCthio.cast <- print(wilcox.test(arsCthio$Count~arsCthio$Classification))

#intI
intI <- subset(x = data.mann, subset = Gene == "intI")
intI.cast <- print(wilcox.test(intI$Count~intI$Classification))

#ClassA
ClassA <- subset(x = data.mann, subset = Gene == "ClassA")
ClassA.cast <- print(wilcox.test(ClassA$Count~ClassA$Classification))

#ClassB
ClassB <- subset(x = data.mann, subset = Gene == "ClassB")
ClassB.cast <- print(wilcox.test(ClassB$Count~ClassB$Classification))

#ClassC
ClassC <- subset(x = data.mann, subset = Gene == "ClassC")
ClassC.cast <- print(wilcox.test(ClassC$Count~ClassC$Classification))

#vanA
vanA <- subset(x = data.mann, subset = Gene == "vanA")
vanA.cast <- print(wilcox.test(vanA$Count~vanA$Classification))

#vanB
vanB <- subset(x = data.mann, subset = Gene == "vanB")
vanB.cast <- print(wilcox.test(vanB$Count~vanB$Classification))

#vanH
vanH <- subset(x = data.mann, subset = Gene == "vanH")
vanH.cast <- print(wilcox.test(vanH$Count~vanH$Classification))

#############################
#EXAMINE CLASS LEVEL CHANGES#
#############################

#order class by phylum
data.annotated.ncbi$class <- factor(data.annotated.ncbi$class, 
                          levels = data.annotated.ncbi$class[order(data.annotated.ncbi$phylum)])

#summarise class level information
data.class <- data.annotated.ncbi %>%
  group_by(Group, Description, Gene, class, Classification, Site, Temp) %>%
  summarise(Class.count = sum(Normalized.Abundance.census))

#prep colors for class diversity
#n <- 64
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#color.class <- print(sample(col_vector, n))
color.class <- c( "#FFD92F", "#F2F2F2", "#CAB2D6", "#CBD5E8", "#F781BF", "#E41A1C", "#FFFF33", "#FDB462", "#FDDAEC", "#CCCCCC", "#FB9A99", "#FFFF99", "#D95F02", "#7570B3", "#CCEBC5", "#8DA0CB", "#FC8D62", "#FFFFB3", "#66C2A5", "#1F78B4", "#E7298A", "#B3CDE3", "#A6761D", "#80B1D3", "#377EB8", "#666666", "#BC80BD", "#FED9A6", "#A6CEE3", "#E78AC3", "#B2DF8A", "#66A61E", "#6A3D9A", "#BEBADA", "#CCEBC5", "#999999", "#D9D9D9", "#E5D8BD", "#FF7F00", "#FDBF6F", "#4DAF4A", "#B15928", "#FDCDAC", "#E5C494", "#666666", "#FFED6F", "#A65628", "#FB8072", "#FFFF99", "#A6D854", "#B3B3B3", "#33A02C", "#8DD3C7", "#F0027F", "#F1E2CC", "#386CB0", "#BF5B17", "#7FC97F", "#F4CAE4", "#FF7F00", "#E31A1C", "#E6AB02", "#FFF2AE", "#BEAED4")

#order genes by group
data.class$Gene <- factor(data.class$Gene, 
                           levels = data.class$Gene[order(data.class$Group)])


#plot arsenic resistance genes with class info
(asrg.gene.bar.class <- ggplot(subset(data.class, Group == "ArsenicResistance"), aes(x = Site, 
                                            y = Class.count*100, fill = class)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = color.class) +
    theme_classic(base_size = 10) +
    ylab("Genome Equivalents with Gene (%)") +
    facet_wrap(~ Gene) +
    theme(axis.text.x = element_text(angle = 90, size = 10, hjust=0.95,vjust=0.2)))

#save plot
ggsave(asrg.gene.bar.class, filename = paste(wd, "/figures/class.abundance.by_AsRG.png", sep=""), width = 10)

#plot arsenic resistance genes with class info
(abrg.gene.bar.class <- ggplot(subset(data.class, Group == "AntibioticResistance"), aes(x = Site, 
                                                                                     y = Class.count*100, fill = class)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = color.class) +
    theme_classic(base_size = 10) +
    ylab("Genome Equivalents with Gene (%)") +
    facet_wrap(~ Gene) +
    theme(axis.text.x = element_text(angle = 90, size = 10, hjust=0.95,vjust=0.2)))

#save plot
ggsave(abrg.gene.bar.class, filename = paste(wd, "/figures/class.abundance.by_abRG.png", sep=""), width = 10)

###########################
#EXAMINE OTU LEVEL CHANGES#
###########################
#remember 1 row in data.annotated.ncbi = 1 OTU (clustered at lowest taxonomy; for most genes this is ~0.1)

(asrg.otu.plot <- ggplot(subset(data.annotated.ncbi, 
                                Group == "ArsenicResistance"),
                    aes(x = Gene, y = Normalized.Abundance.census, 
                        color = phylum, shape = Classification)) +
   geom_jitter(alpha = 0.5) +
   scale_color_manual(values = color) +
   ylab("Genome equivalents with OTU (%)") +
   theme_classic())

#make colors for AB-rg OTUs
col13=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "black", "#AA4499")

#plot ARG OTUs
(abrg.otu.plot <- ggplot(subset(data.annotated.ncbi, 
                                Group == "AntibioticResistance"),
                         aes(x = Gene, y = Normalized.Abundance.census, 
                             color = phylum, shape = Classification)) +
    geom_jitter(alpha = 0.7, width = 0.3, size = 3) +
    scale_color_manual(values = tol12qualitative) +
    ylab("Genome equivalents with OTU (%)") +
    theme_classic())

#extract high abundance organisms (OTUs)
high.abund <- data.annotated.ncbi %>%
  group_by(Description, Group, Gene, Organism) %>%
  summarise(Organism.abundance = length(Organism))

#make list of highest abundance organisms
high.abund.list <- high.abund[which(high.abund$Organism.abundance > 12),]

#make high abundance annotated dataset
data.annotated.ncbi.abundt <- data.annotated.ncbi[which(data.annotated.ncbi$Organism %in% high.abund.list$Organism),] 

#look at acr3 OTUs based on organism name
ggplot(subset(data.annotated.ncbi.abundt, Gene == "acr3"), 
       aes(x = Organism, y = Normalized.Abundance.census, 
                                shape = Classification, 
                                color = Temp)) +
  geom_jitter(alpha = 0.9, size = 3, width = 0.1) +
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                        guide_legend(title="Temperature (°C)")) +
  facet_wrap(~Gene) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust=0.95,vjust=0.2))
  



