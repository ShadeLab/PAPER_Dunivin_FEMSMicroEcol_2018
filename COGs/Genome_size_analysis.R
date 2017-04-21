#load packages
library(reshape2)
library(ggplot2)
library(psych)

#print working directory for future references
wd=print(getwd())

##View trend of average genome size 
#read in microbe census data
data=read.delim(file = paste(wd, "/data/microbe_census.txt", sep=""))

#read in site information
meta=read.delim(file = paste(wd, "/data/Cen_temp_As.txt", sep=""), header = FALSE)
colnames(meta)=c("Site", "Temperature", "As", "Classification")

#join datasets to annotate microbe census data
data=inner_join(data, meta, by="Site")

#plot data
(AGS <- ggplot(data, aes(x=Temperature, y=AGS/1000000)) +
  geom_smooth(method="lm") + 
  geom_point(aes(shape=Classification), size=3) +
  ylab("Average genome size (Mbp)") +
  theme_classic(base_size = 12))

#save plot
ggsave(AGS, filename = paste(wd, "/figures/AGS.temp.png", sep=""))

##calculate correlation
#subset data
cor=cor.test(data$AGS, data$Temperature)
