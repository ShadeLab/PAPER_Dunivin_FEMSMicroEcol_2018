#load required packages
library(tidyverse)
library(psych)
library(reshape2)

#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis
wd <- print(getwd())

#make color pallette for Centralia temperatures
GnYlOrRd <- colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)

#setwd to diversity analysis
setwd(paste(wd, "/diversity_analysis", sep = ""))

#now change wd obj
wd <- print(getwd())

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data/gc_counts", sep = ""))

#read in abundance data
names <- list.files(pattern = "*gc_out.txt")
gc <- do.call(rbind, lapply(names, function(X) {
  data.frame(read_delim(X, delim = "\t"))}))

#read in hmm length data
length <- read_delim("hmm.lengths.txt", delim = "\t")

#move back up a directory to proceed with analysis
setwd("../../")
wd <- print(getwd())

#adjust column names of data file
colnames(gc) <- c("id", "Perc.GC", "Total.Count", "G.Count", 
                  "C.Count", "A.Count", "T.Count")
#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_FULL_map.txt", 
                                    sep=""), sep=" ", header=TRUE))

#capitolize site name
gc$id <- gsub("cen", "Cen", gc$id)

#change arsC_glut and arsC_thio to not include _
gc$id <- gsub("arsC_glut", "arsCglut", gc$id)
gc$id <- gsub("arsC_thio", "arsCthio", gc$id)


#tidy data and spread ids
gc.tidy <- gc %>%
  separate(col = id, into = c("Site", "Gene", "Contig1", "Contig1n", 
                              "Contig2", "Contig2n"), sep = "_") %>%
  unite(Contig, Contig1:Contig2n, remove = TRUE) %>%
  left_join(meta, by = "Site") %>%
  left_join(length, by = "Gene") %>%
  mutate(Length.bp = Length.aa * 3) %>%
  group_by(Site, Gene) 

# temporarily remove cen13
gc.tidy <- gc.tidy[-which(gc.tidy$Site == "Cen13"),]

#remove all rows that are not at least 90% of nucleotide length
gc.tidy.long <- gc.tidy[which(gc.tidy$Total.Count > 0.90*gc.tidy$Length.bp),]

#order sites based on temperature
gc.tidy.long$Site <- factor(gc.tidy.long$Site, 
                                   levels = gc.tidy.long$Site[order(gc.tidy.long$SoilTemperature_to10cm)])

#plot data by per GC
ggplot(gc.tidy.long, aes(x = Site, y = Perc.GC)) +
  geom_boxplot() +
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                        guide_legend(title="Temperature (°C)")) +
  facet_wrap(~Gene) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust=0.95,vjust=0.2))

#plot data by per GC
ggplot(gc.tidy.long, aes(x = Classification, y = Perc.GC)) +
  geom_boxplot() +
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                        guide_legend(title="Temperature (°C)")) +
  facet_wrap(~Gene) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust=0.95,vjust=0.2))

#summarise dataset
gc.summary <- gc.tidy.long %>%
  summarise(Average.perc.gc = mean(Perc.GC), SD = sd(Perc.GC)) %>%
  left_join(meta, by = "Site")

#plot data by temperature
ggplot(gc.summary, aes(x = Classification, y = Average.perc.gc)) +  
  geom_boxplot() + 
  geom_jitter(aes(color = SoilTemperature_to10cm)) +
  scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", 
                        guide_legend(title="Temperature (°C)")) +
  facet_wrap(~Gene) +
  theme_classic()


ggplot(gc.summary, aes(x = SoilTemperature_to10cm, y = Average.perc.gc)) +
  geom_errorbar(aes(ymin=Average.perc.gc-SD, 
                  ymax=Average.perc.gc+SD)) +
  geom_point() +
  facet_wrap(~Gene) +
  theme_classic()

#make wide dataset for correlations
gc.cast <- gc.summary %>%
  dcast(Site ~ Gene, value.var = "Average.perc.gc") %>%
  left_join(meta, by = "Site") %>%
  select(Site:vanZ, SoilTemperature_to10cm)

#select genes with presence in more than 3 sites
gc.cast <- gc.cast %>%
  select(acr3, arsA, arsCglut, arsD, arsM, intI, rplB,
         sul2, vanA, vanH, vanX, vanZ, SoilTemperature_to10cm)


#test gc correlations for each gene
gc.corr <- print(corr.test(gc.cast[,1:13], method = "spearman"), 
                 short = FALSE)
###acr3 GC content is significantly correlated w temp

##############################
#Examine sequence dist length#
##############################
#order sites based on temperature
gc.tidy$Site <- factor(gc.tidy$Site, 
                            levels = gc.tidy$Site[order(gc.tidy$SoilTemperature_to10cm)])

#plot percent total length by site
ggplot(gc.tidy, aes(x = Site, y = Total.Count/Length.bp)) +
  geom_boxplot() +
  geom_jitter(size = 0.5, alpha = 0.5) +
  facet_wrap(~Gene) +
  ylim(0,1.05) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, 
                                   hjust=0.95,vjust=0.2)) 

##make column of sequence identifiers
#uncapitolize site name
gc.tidy.long$Site <- gsub("Cen", "cen", gc.tidy.long$Site)

#change arsC_glut and arsC_thio to include _
gc.tidy.long$Gene <- gsub("arsCglut", "arsC_glut", gc.tidy.long$Gene)
gc.tidy.long$Gene <- gsub("arsCthio", "arsC_thio", gc.tidy.long$Gene)

#make column of contig names
contig.names <- gc.tidy.long %>%
  unite(Contig.Names, Site:Contig, sep = "_", remove = FALSE) %>%
  ungroup() %>%
  select(Contig.Names)
  
#extract list of near-full length sequences
write.table(contig.names, file = paste(wd, "/output/full.length.contigs.txt", 
                                       sep = ""), 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

