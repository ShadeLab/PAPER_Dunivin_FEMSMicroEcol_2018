library(ggplot2)
library(reshape2)
library(dplyr)
library(stats)

#print working directory for future references
wd=print(getwd())

#read in data
data=read.delim(file = paste(wd, "/data/COG_data.txt", sep=""))

#make COG ID's the row names
row.names(data)=data$COGID
data$COGID=NULL

#make a record of functional name and COG ID for later
funct=subset(data, select=Func_name)

#remove functional information
data=subset(data, select=-Func_name)

#read in genome info
size=read.delim(file = paste(wd, "/data/microbe_census.txt", sep=""))

#read in site info
info=read.delim(file = paste(wd, "/data/Cen_temp_As.txt", sep=""), header = FALSE)

#add colnames to site
colnames(info)=c("Site", "Temp", "As", "History")

#tidy cog data
data$COGID=rownames(data)
tidy=melt(data, id.vars="COGID", variable.name="Site", value.name = "Abundance")

#fix sit id naming in cog data
tidy$Site=gsub("C", "Cen", tidy$Site)

#join data with info
joined=inner_join(size, info)
final=inner_join(tidy, joined, by="Site")

#set Abundance to numeric
final$Abundance=as.numeric(final$Abundance)

#normalize abundance
final$norm=final$Abundance/final$GE

#wide dataset
cast=acast(final, Site ~ COGID, id=c("Site", "COGID"), value.var = "norm")
cast2=cbind(cast, joined$Temp)

#correlation
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

##examine single copy genes
#list genes
genes=c("COG0016","COG0048","COG0049","COG0051","COG0052","COG0072","COG0080","COG0081","COG0087","COG0088", "COG0090", "COG0091","COG0092","COG0093","COG0094","COG0096","COG0097","COG0098","COG0099","COG0100","COG0103","COG0127","COG0149","COG0164","COG0184","COG0185","COG0186","COG0197","COG0200","COG0244","COG0256","COG0343","COG0481","COG0504","COG0532","COG0533","COG0541")

#select COGs Tringe lists as single copy
single=final[which(final$COGID %in% genes),]

#group
single=group_by(single, COGID)

#average count
avglength=summarise(single, n=length(norm), average=mean(norm))

#calculate odds ratio
single$odds.ratio=single$norm/avglength$average


#plot odds ratios for each COG (single copy)
(or=ggplot(single, aes(x=COGID, y=odds.ratio)) +
  geom_boxplot() +
  geom_jitter(height = 0, aes(color=Site, shape=History), size=1.5) +
  ylim(0,2) +
  coord_flip())

#save plot
ggsave(or, filename = paste(wd, "/figures/GE.norm.odds.ratios.png", sep=""))

#plot odds ratio v. temperature 
(or.temp=ggplot(single, aes(x=Temp, y=odds.ratio, color=COGID)) +
  geom_point(size=1) +
  stat_smooth(method=lm, alpha=0.01, size=0.25) +
  ylim(0.25, 1.6))

#save plot
ggsave(or.temp, filename = paste(wd, "/figures/GE.norm.odds.ratio.v.temp.png", sep=""))

#check correlations
single.r=slim[which(slim$COGID %in% genes),]











