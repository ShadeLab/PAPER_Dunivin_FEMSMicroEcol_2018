# Assembly assessment
### Taylor Dunivin

### Purpose: 
  * Examine kmer abundance distribution
  * Calculate length (nucl) statistics
  * Count number of sequences
  * Calculate e-values for all sequences (make sure >10^-5)
  
---

## Kmer abundance distribution
We expect curve will have a power law distribution. Deviations from this should be considered further. 
```
#load R
module load GNU/4.9
module load R/3.3.0
R

#load required packages
library(ggplot2)

#read in taxonabundance file
kmer=read.table(list.files(pattern = "_abundance.txt"), header=TRUE)

#plot dist
plot=ggplot(kmer, aes(x=kmer_abundance, y=frequency)) +
  geom_point() +
  labs(x="kmer abundance", y="Frequency")
  
#save plot
ggsave("kmerabundancedist.png", plot=last_plot(), width=4, height=4)
```

## Nucleotide length statistics
```
#remove stats info from framebot file
grep "STATS" *_framebot.txt > framebotstats.txt

#load R
module load GNU/4.9
module load R/3.3.0
R

#load proper modules
library(dplyr)

#read in stats on length
stats=read.table("framebotstats.txt", header=FALSE)

#calculate statistics
results=summarise(stats, ProteinContigClusters.99=length(stats$V4),AverageLength=mean(stats$V4),MedianLength=median(stats$V4), MinLength.bp=min(stats$V4), MaxLength.bp=max(stats$V4), MaxPercentIdentity=max(stats$V6), MinPercentIdentity=min(stats$V6), AveragePercentIdentity=mean(stats$V6))

#save results
write.table(results, "/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/k45/arsC_thio/cluster/stats.txt", row.names=FALSE)
