#load packages
library(reshape2)
library(ggplot2)
library(psych)

##View trend of average genome size 
#read in microbe census data
data=read.delim("microbe_census.txt")

#read in site information
meta=read.delim("Cen_temp_As.txt", header=FALSE)
colnames(meta)=c("Site", "Temperature", "As", "Classification")

#join datasets to annotate microbe census data
data=inner_join(data, meta, by="Site")

#plot data
ggplot(data, aes(x=Temperature, y=AGS/1000000)) +
  geom_point(aes(shape=Classification), size=3) +
  geom_smooth(method="lm") + 
  ylab("Average genome size (Mbp)") +
  theme_classic(base_size = 12)

##calculate correlation
#subset data
cor=cor.test(data$AGS, data$Temperature)
