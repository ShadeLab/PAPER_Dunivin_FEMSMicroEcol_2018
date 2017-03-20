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
module load R/3.3.0
R
library("data.table")
ga=fread("gene2accession", header=TRUE, select=c("tax_id", "GeneID", "protein_accession.version", "protein_gi", "genomic_nucleotide_accession.version", "genomic_nucleotide_gi"))
```
