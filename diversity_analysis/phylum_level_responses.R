#load
library(dplyr)
library(ggplot2)
library(reshape2)

#read in abundance data
data=read.csv(file = paste(wd, "/data/rplB_taxon_abundance_cen.csv", sep=""))

#group data
grouped=group_by(data, Site, Taxon)

#decast for abundance check
dcast=acast(grouped, Taxon ~ Site, value.var = "Fraction.Abundance")

#call na's zeros
dcast[is.na(dcast)] =0

#order based on abundance
order.dcast=dcast[order(rowSums(dcast),decreasing=TRUE),]

#melt data
melt=melt(order.dcast,id.vars=row.names(order.dcast), variable.name= "Site", value.name="Fraction.Abundance" )

#adjust colnames of melt
colnames(melt)=c("Taxon", "Site", "Fraction.Abundance")

#read in metadata
meta=read.delim(file = paste(wd, "/data/Centralia_full_map.txt", sep=""), sep=" ")

#join metadata with regular data
joined=inner_join(melt, meta, by="Site")

#average by fire history
#group data
grouped2=group_by(joined, Taxon, Classification)

#calculate
history=summarise(grouped2, N=length(Fraction.Abundance), Average=mean(Fraction.Abundance))

#plot
phylum.plot=(ggplot(history, aes(x=Taxon, y=Average)) +
  geom_point(size=2) +
  facet_wrap(~Classification, ncol = 1) +
  labs(x="Phylum", y="Mean relative abundance")+
  theme(axis.text.x = element_text(angle = 90, size = 12, hjust=0.95,vjust=0.2))) 

#save plot
ggsave(phylum.plot, filename = paste(wd, "/figures/phylum.responses.png", sep=""), width = 5, height = 5)
