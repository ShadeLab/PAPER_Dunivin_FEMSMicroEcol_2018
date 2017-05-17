#load packages
library(dplyr)
library(reshape2)
library(ggplot2)

#read in data
data=read.delim(file = paste(wd, "/data/COG_data.txt", sep=""))

#make file of single copy gene cogs
genes=c("COG0016","COG0048","COG0049","COG0051","COG0052","COG0072","COG0080","COG0081","COG0087","COG0088", "COG0090", "COG0091","COG0092","COG0093","COG0094","COG0096","COG0097","COG0098","COG0099","COG0100","COG0103","COG0127","COG0149","COG0164","COG0184","COG0185","COG0186","COG0197","COG0200","COG0244","COG0256","COG0343","COG0481","COG0504","COG0532","COG0533","COG0541")
rpoB="COG0085"

#select COGs Tringe lists as single copy
single=data[which(data$COGID %in% genes),]
rpoB <- data[which(data$COGID %in% rpoB),]

#melt the data
single <- melt(single, id.vars = c("COGID", "Func_name"), variable.name = "Site", value.name = "Count")
rpoB <- melt(rpoB, id.vars = c("COGID", "Func_name"), variable.name = "Site", value.name = "Count")
data <- melt(data, id.vars = c("COGID", "Func_name"), variable.name = "Site", value.name = "Count")

#summarise the data
single.sum <- single %>%
  group_by(Site) %>%
  summarise(Average.single = mean(Count))
rpoB.sum <- rpoB %>%
  group_by(Site) %>%
  summarise(Average.rpoB = mean(Count))

##############
#AsRG CAREER##
##############
data.rpoB <- data %>%
  left_join(rpoB.sum, by = "Site") %>%
  mutate(Normalized.count = Count / Average.rpoB)

#select arsB/ACR3 cogs only
arsB <- data.rpoB[which(data$COGID == "COG0798"),]

#read in meta data
meta=data.frame(read.delim(file = "/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/diversity_analysis/data/Centralia_JGI_map.txt", sep=" ", header=TRUE))

#change "cen" to "Cen" so it matches outside data
arsB$Site <- gsub("C", "Cen", arsB$Site)
arsB <- arsB %>%
  left_join(meta, by = "Site") %>%
  select(Site, Normalized.count, Classification)

#read in other soil ref
other <- data.frame(read.delim(file = "/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/COGs/data/other_ref.txt"))

#combine datasets
arsB <- rbind(arsB,other)
#order based on abundance
arsB$Site <- factor(arsB$Site, levels = arsB$Site[order(arsB$Classification)])

arsB.o=arsB[order(arsB$Classification,decreasing=FALSE),]
#plot
(arsB.plot <- ggplot(arsB, aes(x = Site, y = Normalized.count, 
                               fill = Classification)) +
    geom_bar(stat = "identity", color = "black") +
    ylab("Normalized COG count") +
    scale_fill_manual(values = c("firebrick2", "yellow1", "green3", "dodgerblue")) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2)))

ggsave(arsB.plot, filename = paste(wd, "/figures/arsB.plot.eps", sep=""))
#####
#done
#####

#regroup w Site 
reshaped=group_by(single, COGID, Site)

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

