#load required packages
library(tidyverse)
library(reshape2)
library(stringr)

#set up environment
setwd("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/diversity_analysis/")
wd <- paste(getwd())

#read in label data
labels <- read_delim("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/phylogenetic_analysis/tree_data/dissimilatory/dissimilatory_alignment_muscle_0.1_labels_short.txt", delim = "\t", col_names = FALSE)

#remove rows that do not contain "seq_id="
labels <- labels[grep(pattern = "seq_id=", labels$X1),]

#remove > from data
labels$X1 <- gsub(">", "", labels$X1)

#tidy label information
labels.tidy <- labels %>%
  separate(col = X1, into = c("OTU", "Contig"), 
           sep = " seq_id=") %>%
  separate(col = Contig, into = c("Site", "Gene"), 
           sep = "_")

# Add otu to beginning of otu column
#labels.tidy$OTU <- paste("OTU_", labels.tidy$OTU, sep="")

#list otus that are gene matches
#mixed.aioA <- c("OTU_0013", "OTU_0031", "OTU_0034", "OTU_0036", "OTU_0037", "OTU_0042", "OTU_0043", "OTU_0049", "OTU_0051")
#mixed.arrA <- c("OTU_0001", "OTU_0002", "OTU_0003", "OTU_0004", "OTU_0005", "OTU_0006", "OTU_0008")
#mixed.arxA <- c("OTU_0011", "OTU_04")

#remove OTUs that are gene matches
#labels.tidy <- labels.tidy[-which(labels.tidy$Gene == "aioA" &
#                                    labels.tidy$OTU %in% mixed.aioA),]
#labels.tidy <- labels.tidy[-which(labels.tidy$Gene == "arrA" &
#                                    labels.tidy$OTU %in% mixed.arrA),]
#labels.tidy <- labels.tidy[-which(labels.tidy$Gene == "arxA" &
#                                    labels.tidy$OTU %in% mixed.arxA),]

#set working directory different
setwd("../diversity_analysis/data/0.1_clust/")
wd2 <- paste(getwd())

#list filenames of interest
filenames <- c("arrA_rformat_dist_0.1.txt", "arxA_rformat_dist_0.1.txt",
               "aioA_rformat_dist_0.1.txt")

#make dataframes of all OTU tables
for(i in filenames){
  filepath <- file.path(paste(wd2,i,sep="/"))
  assign(gsub("_rformat_dist_0.1.txt", "", i), data.frame(id = basename(i), read.delim(filepath,sep = "\t")))
}

#tidy data before combining
aioA <- melt(aioA, id.vars = c("id", "X"), variable.name = "OTU", value.name = "Abundance")
arrA <- melt(arrA, id.vars = c("id", "X"), variable.name = "OTU", value.name = "Abundance")
arxA <- melt(arxA, id.vars = c("id", "X"), variable.name = "OTU", value.name = "Abundance")

#join data
data <- rbind(aioA, arrA, arxA)

#clean up information
data <- data %>%
  separate(id, into = c("Gene", "junk"), sep = "_") %>%
  rename(Site = X) %>%
  select(-junk)

setwd("../../../diversity_analysis/output/")
wd2 <- paste(getwd())
#read in census data
rplB <- read_delim(file = paste(wd2, "rplB.summary.scg.txt", sep = "/"),
                     delim = " ")

#fix wd
setwd(wd)
#read in metadata
meta <- data.frame(read.delim(paste(wd, "/data/Centralia_FULL_map.txt", 
                                    sep=""), sep=" ", header=TRUE))

#add census data for normalization purposes
table.census <- data %>%
  left_join(rplB, by = "Site") %>%
  mutate(normalized.abundance = Abundance/rplB)

##fix otu labels to match label labels
#remove otu
table.census$OTU <- gsub("OTU_", "", table.census$OTU)
table.census$OTU <- as.numeric(table.census$OTU)

#add leading zero to 4 digits
table.census$OTU <- sprintf("%04d", table.census$OTU)

#replace gene names with numbers
table.census$Gene <- gsub("aioA", "GENE1", table.census$Gene)
table.census$Gene <- gsub("arrA", "GENE2", table.census$Gene)
table.census$Gene <- gsub("arxA", "GENE3", table.census$Gene)

#complete new otu name
table.census$OTU <- paste(table.census$Gene, table.census$OTU, 
                          sep = "")
#spread data
table.normalized <- acast(table.census, OTU~Site, value.var = normalized.abundance)

#rename site column
table.normalized$Site <- table.normalized$`table.census$Site`
table.normalized <- table.normalized[,-1]

#order based on temperature
table.normalized <- table.normalized %>%
  left_join(meta, by = "Site") %>%
  arrange(Site, desc(SoilTemperature_to10cm)) %>%
  select(1:43, 62, 47)

table.normalized$Site <- factor(table.normalized$Site, 
                                levels = table.normalized$Site[order(meta$SoilTemperature_to10cm)])

#transform otu table
table.normalized.t = setNames(data.frame(t(table.normalized[,-c(45,44)])), 
                              table.normalized[,45])
#make OTUs a column
table.normalized.t$OTU <- rownames(table.normalized.t)

#temporarily remove otu_ from otu column
table.normalized.t$OTU <- gsub("OTU_", "", table.normalized.t$OTU)
table.normalized.t$OTU <- as.numeric(table.normalized.t$OTU)

#add leading zero to 4 digits
table.normalized.t$OTU <- sprintf("%04d", table.normalized.t$OTU)

#add OTU to otu label
table.normalized.t <- within(table.normalized.t, Name[Name == 'John Smith' & State == 'WI'] <- 'John Smith1')

#remove weird space after otu number
table.normalized.t$OTU <- gsub(" ", "", table.normalized.t$OTU)

#join abundance with name
label.abund <- left_join(labels.tidy, table.normalized.t, by = "OTU")

#remove rows with NA values
label.abund <- label.abund[!is.na(label.abund$Cen01),]


#remove otu from otu
label.abund$OTU <- gsub("OTU_", "", label.abund$OTU)

label.abund$OTU <- paste(label.abund$Gene, label.abund$OTU, sep = "")

#remove OTU information
label.abund <- select(label.abund, -c(Label,Site,Gene))

#save file
write.csv(label.abund, "/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/phylogenetic_analysis/tree_data/dissimilatory/dissimilatory_abund_label.csv", quote = FALSE, row.names = FALSE)




