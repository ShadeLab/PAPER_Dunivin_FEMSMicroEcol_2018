#load packages
library(dplyr)
library(reshape2)
library(ggplot2)

#read in data
data=read.delim(file = paste(wd, "/data/COG_data.txt", sep=""))

#make file of single copy gene cogs
genes=c("COG0016","COG0048","COG0049","COG0051","COG0052","COG0072","COG0080","COG0081","COG0087","COG0088", "COG0090", "COG0091","COG0092","COG0093","COG0094","COG0096","COG0097","COG0098","COG0099","COG0100","COG0103","COG0127","COG0149","COG0164","COG0184","COG0185","COG0186","COG0197","COG0200","COG0244","COG0256","COG0343","COG0481","COG0504","COG0532","COG0533","COG0541")

#select COGs Tringe lists as single copy
single=data[which(data$COGID %in% genes),]

#read in list of gene lengths
length=read.delim("COGhmmLength.txt", header=FALSE)

#relabel length file based on data
colnames(length)=c("COGID", "HMMlength")

#join length with data
single=inner_join(single, length)

#melt the data
reshaped=melt(single, id.vars = c("COGID", "Func_name", "HMMlength"), variable.name = "Site", value.name = "Count")

#read in gbase data
gb=read.delim("Cen_filesize.txt", header=TRUE)

#add gbases info
reshaped=inner_join(reshaped, gb, by="Site")

#calculate normalized values (by length)
reshaped$LengthNormCount=reshaped$Count/reshaped$HMMlength

#calculate normalized values by gb and length
reshaped$gbNormCount=reshaped$LengthNormCount/reshaped$Gbases

#group the data
reshaped=group_by(reshaped, COGID)

#summarise 
avglength=summarise(reshaped, gbNorm.n=length(gbNormCount), lengthNorm.n=length(LengthNormCount), gbNorm.average=mean(gbNormCount), lengthNorm.average=mean(LengthNormCount))

#regroup w Site 
reshaped=group_by(reshaped, COGID, Site)

#calculate odds ratio
reshaped$length.oddsratio=reshaped$LengthNormCount/avglength$lengthNorm.average
reshaped$gb.oddsratio=reshaped$gbNormCount/avglength$gbNorm.average
  
#plot
ggplot(reshaped, aes(x=COGID, y=length.oddsratio, color=Site)) +
  geom_jitter() +
  coord_flip()

ggplot(reshaped, aes(x=COGID, y=gb.oddsratio, color=Site)) +
  geom_jitter() +
  ylim(0,2) +
  coord_flip()


##plot against
lnc=reshaped[,-c(2,3,5,6,8,9,10)]
gnc=reshaped[,-c(2,3,5,6,7,9,10)]
orig=single[,-c(1,14,15)]

library(tidyr)
lnc.spread=spread(lnc, Site, LengthNormCount)
gnc.spread=spread(gnc, Site, gbNormCount)

library(corrplot)

lnc.spread2=lnc.spread[,-1]

cor=cor(lnc.spread2)
corrplot(cor)

gnc.spread2=gnc.spread[,-1]
gnc.cor=cor(gnc.spread2)
corrplot(gnc.cor)

orig.cor=cor(orig)
corrplot(orig.cor)

#what genes have best odds ratio?
ggplot(reshaped, aes(x=COGID, y=gb.oddsratio)) +
  geom_boxplot() +
  coord_flip()

reshaped3=group_by(reshaped, COGID)
reshaped3summary=summarise(reshaped3, average=mean(gb.oddsratio), n=length(gb.oddsratio), sd=sd(gb.oddsratio))
ggplot(reshaped3summary, aes(x=COGID, y=average)) +
  geom_bar(stat="identity") +
  geom_errorbar(ymin=reshaped3summary$average-reshaped3summary$sd, ymax=reshaped3summary$average+ reshaped3summary$sd) +
  ylim(0,2) +
  coord_flip()

