# TetQ: diverse gene references
### Taylor Dunivin

## Table of contents
* [March 6, 2017](https://github.com/ShadeLab/Xander_arsenic/blob/master/reference_sequences/tetQ_references.md#march-6-2017)
* [March 14, 2017](https://github.com/ShadeLab/Xander_arsenic/blob/master/reference_sequences/tetQ_references.md#march-14-2017)
* [March 15, 2017](https://github.com/ShadeLab/Xander_arsenic/blob/master/reference_sequences/tetQ_references.md#march-15-2017)

## March 6, 2017
### Goals: 
* Obtain sets of relatively high quality sequences
* Dereplicate sequences
* Make tree to detemine diversity of sequences
* Assess (and repeat)

### Gather sequences
* TetQ FunGene database already exists
 * authored by Carlos Rodriguez-Minguela
 * 120695 sequences exist
 
| Protein | Test # | min aa | min HMM coverage (%) | min HMM score |
| --------- | ----- | ---------- | --------- | :-----: |
| TetQ | 1 | 620 | 80 | 900 |

* Results in 233 sequences

### Dereplicate sequences
```
java -Xmx2g -jar /mnt/research/rdp/public/RDPTools/Clustering.jar derep -o derep.tetQ.v1.fa all.seqs.tetQ.v1.ids all.seqs.tetQ.v1.samples fungene_8.8_tetQ_233_unaligned_protein_seqs.fa
```
* 56 unique sequences

### Align sequences
Since there are not many sequences, I did not need to submit a job.
```
muscle -in derep.tetQ.v1.fa -out aligned.derep.tetQ.v1.fa
```
* Longest sequence (aa): 670
* Avg length (aa): 641

### Make tree
Since there are not many sequences, I did not need to submit a job.
```
module load raxml/8.0.6
raxmlHPC -m PROTGAMMAJTTF -p 12345 -s aligned.derep.tetQ.v1.fa -n tetq.v1
```

This gives the following outputs:
* ```RAxML_bestTree.tetq.v1```
* ```RAxML_log.tetq.v1```
* ```RAxML_parsimonyTree.tetq.v1```
* ```RAxML_result.tetq.v1```

I used iTol to visualize the RAxML tree (http://itol.embl.de/)
* bestTree
* ![tetQ bestTree](https://github.com/ShadeLab/Xander_arsenic/blob/master/images/RAxML.best.tetQ.pdf)
* Based on this tree, I will remove CAD55718 as it is very distant from other sequences
* parsimony
 * ![tetQ parismony](https://github.com/ShadeLab/Xander_arsenic/blob/master/images/RAxML.parsimony.tetq.png)
 * This tree also suggests I remove CAD55718
 * I am also wondering whether to remove the cluster of sequences right after CAD55718
  * Only a couple of them are called tetQ in NCBI
  * They are quite different from the other sequences in the tree
  * I will leave them in during the initial analysis but will consider removing them if they appear to cause problems
  
### Make nucleotide sequence file to match derep amino acids
Now that I have reduced the total number of sequences for proteins, I need to do the same for nucleotide sequences of tetQ. This cannot be done using derep because this will result in a different number of sequences. There are a couple of options to do this (considering headings have accession number and descriptor): 
1. Download all nucleotide sequences from FunGene with the above specified cutoffs, find a way to delete all sequences that are not associated with remaining protein sequences. This would involve ignoring the accession number (different for protein and nucleotides) and looking only at the descriptor (the same for nucleotide and amino acid sequences). 
2. Blast remaining amino acid sequences and obtain corresponding nucleotide sequences. 

I am not sure how to accomplish either of these with .fa files. I will troubleshoot this before taking detailed notes. 

## March 14, 2017
### Goals
* Re download and derep TetQ files
* Derep nucleotide files based on dereplicated protein

### Download TetQ sequences
* Use nucleotide GI as sequence identifier for BOTH protein and nucleotide data
* Use same cutoffs as before 
| Protein | Test # | min aa | min HMM coverage (%) | min HMM score |
| --------- | ----- | ---------- | --------- | :-----: |
| TetQ | 1 | 620 | 80 | 900 |
* Results in 235 sequences
* Called 
  * fungene_8.8_tetQ_235_unaligned_protein_seqs_v2.fa
  * fungene_8.8_tetQ_235_unaligned_nucleotide_seqs_v2.fa

### Derep protein sequences 
```
java -Xmx2g -jar /mnt/research/rdp/public/RDPTools/Clustering.jar derep -o derep.tetQ.v2.fa all.seqs.tetQ.v2.ids all.seqs.tetQ.v2.samples fungene_8.8_tetQ_233_unaligned_protein_seqs_v2.fa
```
* I get the following error
```
Processing fungene_8.8_tetQ_235_unaligned_nucleotide_seqs_v2.fa
usage: Dereplicator [options] <id-mapping-out> <sample-mapping-out>
                    <seq-file>[,<qual-file>] ...
 -a,--aligned            Dereplicate aligned sequences
 -f,--formatted          Dereplicate formated (uppercase/- = comparable,
                         lowercase/. = non-comparable) aligned sequences
 -g,--keep-common-gaps   Don't remove common gaps in output sequences
 -m,--model-only <arg>   Dereplicate aligned sequences using mask sequence
 -o,--out <arg>          Write sequences to this file
 -q,--qual-out <arg>     Write quality sequences to this file
 -s,--sorted             Sort sequence by number of members represented
 -u,--unaligned          Dereplicate unaligned sequences
java.lang.IllegalArgumentException: Attempting to add duplicate ids [893723728]
	at edu.msu.cme.pyro.derep.IdMapping.addIds(IdMapping.java:76)
	at edu.msu.cme.pyro.derep.Dereplicator.main(Dereplicator.java:295)
	at edu.msu.cme.pyro.cluster.ClusterMain.main(ClusterMain.java:361)
Error: Attempting to add duplicate ids [893723728]
```
* In TextWrangler, I deleted sequences that had identical GIs (an unfortunate caveat of using GI over accno)
```
java -Xmx2g -jar /mnt/research/rdp/public/RDPTools/Clustering.jar derep -o derep.tetQ.v2.fa all.seqs.tetQ.v2.ids all.seqs.tetQ.v2.samples fungene_8.8_tetQ_233_unaligned_protein_seqs_v2.fa
```
* Results
```
Total sequences: 232
Unique sequences: 55
Dereplication complete: 134
```
* The GI duplication did not have a large effect on final results (55 vs. 56 sequences)
* I will move forward with this method

### Make list of derep GI 
* Output file ```all.seqs.tetQ.v2.ids``` contains GI of all sequences compiled together (replicates listed)
* I want only the first "column" (GI) of each row to derep nucleotide sequences
 1. Remove lines (http://www.miniwebtool.com/remove-line-numbers/) and save as csv file
 2. Extract one column using R
```
module load R
R
#set working directory
setwd("/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/tetQ")
#read in file
a=read.table("all.seqs.tetQ.v2.ids.csv", header=FALSE, sep=",", fill=TRUE, col.names = paste0("all.seqs.tetQ.v2.ids.csv",seq_len(max(count.fields("all.seqs.tetQ.v2.ids.csv", sep = ',')))))

#only keep first column (representative sequence)
b=a$all.seqs.tetQ.v2.ids.csv1

#save file
write.table(b, "first.seqs.tetQ.v2.ids", append= FALSE, quote=FALSE, row.names=FALSE, col.names = FALSE)
 ```

### Derep nucleotide sequences based on protein
I already have input and filter name files. 
```
module load bbmap
filterbyname.sh in=fungene_8.8_tetQ_235_unaligned_nucleotide_seqs_v2.fa out=derep.nucl.tetQ.v2.fa names=first.seqs.tetQ.v2.ids include=t
```
Results
```
Input is being processed as unpaired
Time:               0.263 seconds.
Reads Processed:    234 	0.89k reads/sec
Bases Processed:    450375 	1.71m bases/sec
Reads Out:          55
Bases Out:          105780
```

The remaining number are the number of expected sequences. Woo!
__I now have the appropriate sequences to run Xander on TetQ__

I will make a new folder ```xander.seqs``` to hold copies of the final sequences in xander format (with framebot.fa for protein sequences and nucl.fa for nucleotide sequences. 

## March 15, 2017
### Goal: Prepare gene references/ gene db to run Xander

* Set up working directory
	* Need to edit prepare_gene_ref.sh to fit my requirements
```
#!/bin/bash -login

if [ $# -ne 1 ]; then
        echo Usage: gene 
        exit 1
fi

## THIS MUST BE MODIFIED TO YOUR FILE SYSTEM
JAR_DIR=/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/
REF_DIR=/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Xander_assembler

## NOTE you need to used the modified hmmer-3.0_xanderpatch to build the specialized forward and reverse HMMs for Xander 
hmmer_xanderpatch=/mnt/research/rdp/public/thirdParty/hmmer-3.0_xanderpatch/

gene=$1

## REQUIRED FILES 
# This directory must exists:  RDPTools/Xander_assembler/gene/originaldata 
# In the originaldata directory, you need these files (from FunGene site with option Minimum HMM Coverage at least 80 (%)):
# gene.seeds: small set of protein sequences in FASTA format, used to build forward and reverse HMMs
# gene.hmm: this is the HMM built from gene.seeds using original HMMER3. This will be used to build for_enone.hmm and align contigs after assembly 
# nucl.fa: a large near full length known set used by UCHIME chimera check step
# framebot.fa: a large near full length known protein set to be used to create starting kmers and FrameBot 

## OUTPUTS
## This script will create three files to dir RDPTools/Xander_assembler/gene/:
## for_enone.hmm and rev_enone.hmm that will be used by Xander search step.
## ref_aligned.faa will be used by Xander find starting kmer step


cd ${REF_DIR}/gene_resource/${gene}/originaldata || { echo " directory not found" ;  exit 1; }

## create forward and reverse hmms for Xander.
${hmmer_xanderpatch}/src/hmmalign --allcol -o ${gene}_seeds_aligned.stk ${gene}.hmm ${gene}.seeds
${hmmer_xanderpatch}/src/hmmbuild --enone ../for_enone.hmm ${gene}_seeds_aligned.stk

java -jar ${JAR_DIR}/ReadSeq.jar to-fasta ${gene}_seeds_aligned.stk > ${gene}_seeds_aligned.fasta

python ${JAR_DIR}/Xander_assembler/pythonscripts/reverse.py ${gene}_seeds_aligned.fasta
java -jar ${JAR_DIR}/ReadSeq.jar to-stk -r rev_${gene}_seeds_aligned.fasta rev_${gene}_seeds_aligned.stk

${hmmer_xanderpatch}/src/hmmbuild --enone ../rev_enone.hmm rev_${gene}_seeds_aligned.stk

${hmmer_xanderpatch}/src/hmmalign --allcol -o ref_aligned.stk ../for_enone.hmm framebot.fa
java -jar ${JAR_DIR}/ReadSeq.jar to-fasta ref_aligned.stk > ../ref_aligned.faa

rm *stk ${gene}_seeds_aligned.fasta rev_${gene}_seeds_aligned.fasta
```

* Prepare gene data
	* Add all files to folder ```originaldata```
	* Download and add ```tetQ.hmm``` and ```tetQ.seeds``` by selecting "select seed sequences" on the tetQ fungene page (http://fungene.cme.msu.edu/hmm_detail.spr?hmm_id=33) and then selecting "download hmm" and "download seeds."

* Prepare gene reference
``` ./prepare_gene_ref.sh tetQ ```

* Output files: ```for_enone.hmm```, ```ref_aligned.faa```, ```rev_enone.hmm```

__With this, tetQ references for Xander (V1) are complete__
