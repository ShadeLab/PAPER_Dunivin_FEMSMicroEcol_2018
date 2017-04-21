#load packages
library(dplyr)
library(reshape2)
library(ggplot2)

#print working directory for future references
wd=print(getwd())

#read in data
data=read.delim(file = paste(wd, "/data/COG_data.txt", sep=""))

#read in coverage (size in Gb)
size=read.delim(file = paste(wd, "/data/Cen_filesize.txt", sep = ""))

#read in metadata 
meta=read.delim(file = paste(wd, "/data/Cen_temp_As.txt", sep = ""), header = FALSE)
colnames(meta) <- c("Site", "Temp", "As", "Classification")
meta$Site <- as.character(meta$Site)
meta$Site <- gsub("Cen", "C", meta$Site)

#melt the data
reshaped=melt(data, id.vars = c("COGID", "Func_name"), variable.name = "Site", value.name = "Count")

#add Gbases to data
reshaped <- inner_join(reshaped, size, by = "Site") %>%
  inner_join(., meta, by = "Site")

#list genes used by microbe census
mc.genes <- c("COG0052", "COG0081", "COG0532", "COG0091", 'COG0088', 'COG0090', "COG0103", 'COG0087', "COG0072", "COG0093", "COG0098", "COG0185", "COG0049", "COG0197", "COG0099", "COG0016", "COG0200", "COG0097", "COG0080", "COG0094", "COG0048", "COG0092", "COG0100", "COG0244", "COG0096", "COG0256", "COG0184", "COG0186", "COG0102", "COG0198")

#select COGs from microbeCensus 
single=reshaped[which(reshaped$COGID %in% mc.genes),]

#calculate averages of single copy genes
avg <- single %>%
  group_by(Site) %>%
  summarise(N=length(Count), Average=mean(Count), Median=median(Count), SD=sd(Count))

#join datasets
reshaped <- inner_join(reshaped, avg, by = "Site")

#normalize COGs to single copy genes
reshaped$norm=reshaped$Count/reshaped$Average

#calculate normalized values by gb and length
reshaped$norm.size=reshaped$norm/reshaped$Gbases

##Correlations
#widen dataset
cast=acast(final, Site ~ COGID, id=c("Site", "COGID"), value.var = "norm")

#calculate correlation
corr=cor(joined$Temp, cast, method = "pearson")

#tidy correlated data
slim=melt(corr,value.name = "R")
slim=slim[,-1]
colnames(slim)=c("COGID", "R")

#annotate correlated data
funct$COGID=rownames(funct)
slim=inner_join(slim, funct, by="COGID")

#subset based on correl
neg=slim[which(slim$R<0),]
mpos=subset(slim, slim$R>0 & slim$R<0.8)
vpos=slim[which(slim$R>0.9),]

#subset initial data based on correl
neg.full=final[which(final$COGID %in% neg$COGID),]
mpos.full=final[which(final$COGID %in% mpos$COGID),]
vpos.full=final[which(final$COGID %in% vpos$COGID),]


