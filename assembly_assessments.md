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
I want to double check that all e-values are less than 10^-5. Otherwise, I will increase min length xander parameter. Additionally, the top match should be consistent with the gene of interest
```
#load blast
module load GNU/4.4.5
module load Gblastn/2.28

#make database from diverse gene sequences
makeblastdb -in /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Xander_assembler/gene_resource/GENE/originaldata/nucl.fa  -dbtype nucl -out database

##blast results against db
#tabular format, show seq id, query id (and description), e-value, only show 1 match,
blastn -db database -query NAME_GENE_KMER_final_nucl.fasta -out blast.txt -outfmt "6 qseqid salltitles evalue" -max_target_seqs 1
```

__Output:__ database.nhr, database.nin, database.nsq (database results), blast.txt (blast results)

### Number of reads covering kmers
I want information on the number of reads covering kmers since some studies use this as part of normalization
Note: this is done in a directory above clustering 

```
#load R
module load GNU/4.9
module load R/3.3.0
R

#load packages
library(dplyr)

#read in file 
data=read.table("matchreadlist.txt", header=FALSE)

#count unique in query_id column
reads=summarize(data, UniqueReads=length(unique(data$V1)),TotalReads=length(data$V1))

#write results
write.table(reads, "/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/k45/arsC_thio/cluster/readssummary.txt", row.names=FALSE)
```

### Automation
```
#load blast
module load GNU/4.4.5
module load Gblastn/2.28

#make database from diverse gene sequences
#must input GENE name!!!
makeblastdb -in /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Xander_assembler/gene_resource/GENE/originaldata/nucl.fa  -dbtype nucl -out database

##blast xander results against db
#tabular format, show seq id, query id (and description), e-value, only show 1 match,
blastn -db database -query *_final_nucl.fasta -out blast.txt -outfmt "6 qseqid salltitles evalue" -max_target_seqs 1

#make a list of reads from *match_reads.fa
grep "^>" *match_reads.fa | sed '0~1s/^.\{1\}//g' >matchreadlist.txt

##Extract info for stats calculations in R
#take the framebot stats
grep "STATS" *_framebot.txt > framebotstats.txt

#take the list of matching reads (from sequencing)
grep "^>" *match_reads.fa | sed '0~1s/^.\{1\}//g' >matchreadlist.txt

#start in cluster directory from xander output!
#load R
module load GNU/4.9
module load R/3.3.0
R

#load required packages
library(ggplot2)
library(dplyr)

##COUNT NUMBER OF READ MATCHES
#load packages
library(dplyr)

#read in file 
data=read.table("matchreadlist.txt", header=FALSE)

#count unique in query_id column
reads=summarize(data, UniqueReads=length(unique(data$V1)),TotalReads=length(data$V1))

#write results
write.table(reads, "readssummary.txt", row.names=FALSE)

##KMER ABUND DISTRIBUTION
#read in kmer abund file
kmer=read.table(list.files(pattern = "_abundance.txt"), header=TRUE)

#plot dist
plot=ggplot(kmer, aes(x=kmer_abundance, y=frequency)) +
  geom_point() +
  labs(x="kmer abundance", y="Frequency")
  
#save plot
ggsave("kmerabundancedist.png", plot=last_plot(), width=4, height=4)

##NUCL STATS 
#read in stats on length
stats=read.table("framebotstats.txt", header=FALSE)

#calculate statistics
results=summarise(stats, ProteinContigClusters.99=length(stats$V4),AverageLength=mean(stats$V4),MedianLength=median(stats$V4), MinLength.bp=min(stats$V4), MaxLength.bp=max(stats$V4), MaxPercentIdentity=max(stats$V6), MinPercentIdentity=min(stats$V6), AveragePercentIdentity=mean(stats$V6))

#save results
write.table(results, "stats.txt", row.names=FALSE)

##BLAST STATS
#read in blast results from above
blast=read.delim("blast.txt", header=FALSE)

#write column names based on blast search
colnames(blast)=c("contig", "match", "eval")

#calculate number of low quality sequences along with
#the min, max, mean, and median e values
evalues=summarize(blast, lowq=length(blast[,which(blast$eval>1e-2)]), min=min(blast$eval), max=max(blast$eval), avg=mean(blast$eval), median=median(blast$eval))

#save results
write.table(evalues, "e.values.txt", row.names=FALSE)

```
