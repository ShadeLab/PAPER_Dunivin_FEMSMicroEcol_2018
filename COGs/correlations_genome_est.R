library(ggplot2)
library(reshape2)
library(dplyr)
library(stats)
library(tidyr)

##############
#PREPARE DATA#
##############
#print working directory for future references
wd=print(getwd())

#read in data
data=read.delim(file = paste(wd, "/data/abundance_cog_both_est.tab.xls", sep=""))

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

#join data with info
final=size %>%
  inner_join(info) %>%
  inner_join(.,tidy, by="Site")

#set Abundance to numeric
final$Abundance=as.numeric(final$Abundance)

################
#NORMALIZE DATA#
################
#normalize abundance genome equivalents
final$norm=final$Abundance/final$GE

##normalize abundance to total COG counts
#calculate toal COG count per site
total <- final %>%
  group_by(Site) %>%
  summarise(N=length(Abundance), 
            total=sum(Abundance))

#add total data to final df
final <- final %>%
  left_join(total, by = "Site") 

#normalize abundance to total 
final$anorm=final$Abundance/final$total

#########################################
#CALCULATE CORRELATIONS WITH TEMPERATURE#
#########################################
#widen GE-normalized dataset
norm=acast(final, Site ~ COGID, id=c("Site", "COGID"), value.var = "norm")

#widen total-normalized dataset
anorm=acast(final, Site ~ COGID, id=c("Site", "COGID"), value.var = "anorm")

#correlation based on GE-normalized dataset
corr=cor(info$Temp, norm, method = "pearson")
temp <- matrix(rep(info$Temp, 4631), c(12, 4631))
temp <- info$Temp
norm2 <- t(norm)

coretest.out=NULL
for(i in 1:nrow(norm)){
  results=cor.test(norm[i,],temp)
}


coretest.out=cbind(coretest.out,c(row.names(norm)[i],results$estimate,results$p.value))

#correlation based on total-normalized dataset
acorr=cor(info$Temp, anorm, method = "pearson")

#tidy GE-correlated data
tcorr=melt(corr,value.name = "R")
tcorr=tcorr[,-1]
colnames(tcorr)=c("COGID", "GER")

#tidy total-correlated data
tacorr=melt(acorr,value.name = "R")
tacorr=tacorr[,-1]
colnames(tacorr)=c("COGID", "totR")

#join correl data
corr=inner_join(tcorr, tacorr, by = "COGID")


###########################
#EXAMINE SINGLE COPY GENES#
###########################
#list genes
genes=c("COG0016","COG0048","COG0049","COG0051","COG0052","COG0072","COG0080","COG0081","COG0087","COG0088", "COG0090", "COG0091","COG0092","COG0093","COG0094","COG0096","COG0097","COG0098","COG0099","COG0100","COG0103","COG0127","COG0149","COG0164","COG0184","COG0185","COG0186","COG0197","COG0200","COG0244","COG0256","COG0343","COG0481","COG0504","COG0532","COG0533","COG0541")

#list genes used by microbe census
mc.genes <- c("COG0052", "COG0081", "COG0532", "COG0091", 'COG0088', 'COG0090', "COG0103", 'COG0087', "COG0072", "COG0093", "COG0098", "COG0185", "COG0049", "COG0197", "COG0099", "COG0016", "COG0200", "COG0097", "COG0080", "COG0094", "COG0048", "COG0092", "COG0100", "COG0244", "COG0096", "COG0256", "COG0184", "COG0186", "COG0102", "COG0198")

#select COGs Tringe lists as single copy
single=final[which(final$COGID %in% genes),]

#group
single=group_by(single, COGID)

#average count of GE-normalized cogs
avglength=summarise(single, n=length(norm), average=mean(norm))

#average count of total-normalized cogs
aavglength=summarise(single, n=length(anorm), average=mean(anorm))

#calculate odds ratio GE-normalized cogs
single$odds.ratio=single$norm/avglength$average

#calculate odds ratio total-normalized cogs
single$a.odds.ratio=single$anorm/aavglength$average

#plot odds ratios for each COG (single copy) (GE-normalized)
(or=ggplot(single, aes(x=COGID, y=odds.ratio)) +
  geom_boxplot() +
  geom_jitter(height = 0, aes(color=Site, shape=History), size=1.5) +
  ylim(0,2) +
  xlab("Tringe paper COGS") +
  coord_flip())

#save plot
ggsave(or, filename = paste(wd, "/figures/GE.norm.odds.ratios.png", sep=""))

#plot odds ratios for each COG (single copy) (total-normalized)
(or=ggplot(single, aes(x=COGID, y=a.odds.ratio)) +
    geom_boxplot() +
    geom_jitter(height = 0, aes(color=Site, shape=History), size=1.5) +
    ylim(0,2) +
    xlab("Tringe paper COGS") +
    coord_flip())

#plot odds ratio v. temperature (GE-normalized)
(or.temp=ggplot(single, aes(x=Temp, y=odds.ratio, color=COGID)) +
  geom_point(size=1) +
  stat_smooth(method=lm, alpha=0.01, size=0.25) +
  ylim(0.25, 1.6))

#save plot
ggsave(or.temp, filename = paste(wd, "/figures/GE.norm.odds.ratio.v.temp.png", sep=""))

#plot odds ratio v. temperature (total-normalized)
(or.temp=ggplot(single, aes(x=Temp, y=a.odds.ratio, color=COGID)) +
    geom_point(size=1) +
    stat_smooth(method=lm, alpha=0.01, size=0.25) +
    ylim(0.25, 1.6))

#check correlations
single.r=slim[which(slim$COGID %in% genes),]


#select COGs Tringe lists as single copy
census=final[which(final$COGID %in% mc.genes),]

#group
census=group_by(census, COGID)

#average count
mc.avglength=summarise(census, n=length(norm), average=mean(norm))

#calculate odds ratio
census$odds.ratio=census$norm/mc.avglength$average

#plot odds ratios for each COG (single copy)
(mc.or=ggplot(census, aes(x=COGID, y=odds.ratio)) +
    geom_boxplot() +
    geom_jitter(height = 0, aes(color=Site, shape=History), size=1.5) +
    ylim(0,2) +
    xlab("Microbe census COGS") +
    coord_flip())

#save plot
ggsave(mc.or, filename = paste(wd, "/figures/mc.GE.norm.odds.ratios.png", sep=""))

#plot odds ratio v. temperature 
(mc.or.temp=ggplot(census, aes(x=Temp, y=odds.ratio, color=COGID)) +
    geom_point(size=1) +
    stat_smooth(method=lm, alpha=0.01, size=0.25) +
    ylim(0.25, 1.6))

#save plot
ggsave(mc.or.temp, filename = paste(wd, "/figures/mc.GE.norm.odds.ratio.v.temp.png", sep=""))

#######################
#FUNCTIONAL ASSESSMENT#
#######################

#read in COGIDs with functional letter
cogid=read.delim(file = paste(wd, "/data/cognames2003-2014.txt", sep=""))

#add column names
colnames(cogid)=c("COGID", "func", "name")

#tidy COGs with multiple functional groups
cogid= cogid %>%
  separate(func, c("func", "func2"), sep=1) %>%
  separate(func2, c("func2", "func3"), sep=1) %>%
  melt(id.vars = c("COGID", "name"), measure.vars = c("func", "func2", "func3")) 

#remove unnecessary function column
cogid <- cogid[-3]

#remove rows with empty values
cogid <- cogid[-which(cogid$value == ""),]

#change column names of cogid
colnames(cogid) <- c("COGID", "name", "func")

#join cogid data with data
final <- left_join(final, cogid, by = "COGID")

#read in functional letter key
funct.key=read.delim(file = paste(wd, "/data/fun2003-2014.txt", sep=""))
colnames(funct.key)=c("func", "func.name")

#annotate data with COGID functional groups
final <- final %>%
  left_join(funct.key, by = "func")
  
#group appropriately
final <- group_by(final, func, Site)

#plot
(totplot <- ggplot(final, aes(x = func.name, y = anorm, color = History)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip())

ggplot(final, aes(x = func.name, y = anorm, color = History)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip()
  
(geplot <- ggplot(final, aes(x = func.name, y = norm, fill = History)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip())

###########################
#EXAMINE %GENE VS. %GENOME#
###########################
#annotated temperature correlation data
corr.annotated <- corr %>%
  inner_join(., cogid, by = "COGID") %>%
  inner_join(., funct.key, by = "func")

#group genes based on %gene in total and % gene in genomes
pp <- corr.annotated[which(corr.annotated$GER > 0.6 & corr.annotated$totR > 0.6),]
pp$designation <- "pp"
pn <- corr.annotated[which(corr.annotated$GER < 0.6 & corr.annotated$GER > c(-0.6) & corr.annotated$totR > 0.6),]
pn$designation <- "pn"

pN <- corr.annotated[which(corr.annotated$GER < c(-0.6) & corr.annotated$totR > 0.6),]
pN$designation <- "pN"

nN <- corr.annotated[which(corr.annotated$GER < c(-0.6) & corr.annotated$totR < 0.6 & corr.annotated$totR > c(-0.6)),]
nN$designation <- "nN"

NN <- corr.annotated[which(corr.annotated$GER < c(-0.6) & corr.annotated$totR < c(-0.6)),]
NN$designation <- "NN"

nn <- corr.annotated[which(corr.annotated$GER > c(-0.6) & corr.annotated$GER < 0.6 & corr.annotated$totR < 0.6 & corr.annotated$totR > c(-0.6)),]
nn$designation <- "nn"


desi <- rbind(pp, pn, pN, nN, NN, nn)
sumi <- desi %>%
  group_by(func.name, designation) %>%
  summarise(N=length(designation))

write.csv(sumi, "designations.csv", row.names = FALSE)







