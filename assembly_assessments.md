# Assembly assessment
### Taylor Dunivin

### Purpose: 
  * Examine kmer abundance distribution
  * Calculate length (nucl) statistics, % identity statistics, and # at 99%
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
__Output:__ kmerabundancedist.png

## Statistics
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
```
__Output:__ stats.txt

## BLAST e-values
I want to double check that all e-values are less than 10^-5. Otherwise, I will increase min length xander parameter
```
#load blast
module load GNU/4.4.5
module load Gblastn/2.28

#make database from diverse gene sequences
makeblastdb -in /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Xander_assembler/gene_resource/arsC_thio/originaldata/nucl.fa  -dbtype nucl -out database

##blast results against db
#tabular format, show seq id, query id, e-value, only show 1 match,
blastn -db database -query cen10_arsC_thio_45_final_nucl.fasta -out blast.txt -outfmt "6 qseqid sseqid evalue" -max_target_seqs 1
```
__Output:__ database.nhr, database.nin, database.nsq (database results), blast.txt (blast results)
