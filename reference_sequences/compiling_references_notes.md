# Compiling gene references (general protocol notes)
### Taylor Dunivin


## March 20, 2017

### Description
* I have been trying to compile gene references for multiple genes at once
* This is not going well since I keep running into issues with cooridinating protein and nucleotide replication
* I would like to have a single script to apply to all genes 
* This page will detail my notes on how to go about this

### Current issue
I need protein and corresponding nucleotide sequences for each gene of interest
 * I will only derep nucleotide sequences (will set up script so that either nucl or aa derep can be accomplished)
 * I need to derep the other set of sequences BASED ON this derep
  * Nucl and aa sequences have _different_ accession numbers
  * GI numbers will not work since one GI (organism) can have multiple copies of a gene
  * I also would like taxaIDs so that when I make a tree (final analysis), I will be able to assess the diversity of my reference sequences

### Use "gene2accession" file from NCBI 
This file shows nucleotide accession, protein accession, and taxa_id information for all of NCBI
* I have one concern that it is not updated enough, but I will cross that bridge if it comes
* This file is very large >5GB but I want to use it in R to extrapolate my accnos of interest with semi_join
* I will only read columns of interest to reduce the file size

```
module load GNU/4.9
module load R/3.3.0
R

#load data.table since it is good for reading in large datasets
library("data.table")

##read in gene2accno data
#read in file with accno information for protein and nucleotide sequences; only read necessary columns (save time/memory) (ONCE; ~10min)
ga=fread("gene2accession", header=TRUE, select=c("#tax_id", "GeneID", "protein_accession.version", "protein_gi", "genomic_nucleotide_accession.version", "genomic_nucleotide_gi"))

#remove version numbers from ga file (FunGene does not give .version in its accnos)
ga$protein_accession.version=gsub("\\.[0-2]$","",ga$protein_accession.version)
ga$genomic_nucleotide_accession.version=gsub("\\.[0-2]$","",ga$genomic_nucleotide_accession.version)

#since it is a time consuming process, write ga to a real file so that you only have to make it from the full file (10 extra columns) once
write.table(ga, file="/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/intI/ga.txt", col.names=T, row.names=F)
```

### Need to get sequence information in appropriate format
Do not care about sequences in this step
Want to get sequence identifier only AND separate accession number from descriptions (but keep the description)
{insert steps to get seq names}
```
#remove sequence identifiers by extracting all lines starting with >
grep '^>' fungene_8.9_intI_9418_unaligned_nucleotide_seqs_v1.fa > nucl.seq.id.v1.txt

#remove > (not part of identifier)
sed '0~1s/^.\{1\}//g' nucl.seq.id.v1.txt >>nucl.seq.id.v1.final.txt
```

Now that I have a file with sequence identifier info only, I will load it into R and extract gene accession info based on it 

```
module load GNU/4.9
module load R/3.3.0
R

##prep nucleotide sequence names
#read in sequence identifiers from nucl fasta (from steps above) 
nucl=read.csv("nucl.seq.id.v1.final.txt", col.names = paste0("V",seq_len(8)), fill = TRUE)

#remove location information in nucl (dont need it) so that only accno is in col 1
nucl$V1=gsub(" location=.*", "", nucl$V1)

#name columns in nucl to correspond with ga data
colnames(nucl)=c("genomic_nucleotide_accession.version", "organism", "definition", "definition2", "definition3", "definition4")

##extract relevant ga info based on nucl
library(dplyr)
ga.nucl=semi_join(ga, nucl, by=genomic_nucleotide_accession.version)
```
