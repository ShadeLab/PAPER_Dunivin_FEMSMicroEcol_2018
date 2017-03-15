# gyrB diverse gene references
### Taylor Dunivin
## March 6, 2017
### Goals: 
* Obtain a set of relatively high quality GyrB sequences
* Dereplicate sequences
* Make tree to detemine diversity of sequences
* Assess (and repeat)

## Downloading from FunGene
* GyrB FunGene database already exists
  * authored by Zarraz May-Ping Lee
  * 138459 existing sequences
* Need to narrow down based on
  1. HMM coverage - typically don't go below 80%
  2. Length (aa) - GyrB E. coli is 804aa
  3. Score - based on where there is a large score drop off

| Protein | Test # | min aa | min HMM coverage (%) | min HMM score |
| --------- | ----- | ---------- | --------- | :-----: |
| GyrB | 1 | 800 | 90 | NA |

* The above talbe shows the filters I'm currently trying. I did not pick a minimum score since there was not a clear drop off. 

## Dereplicate sequences
* Use the RDP's derep function fom Clustering
* Example code

```
java -Xmx2g -jar /path/to/Clustering.jar derep -o derep.gyrB.fa all.seqs.gyrB.ids all.seqs.gyrB.samples file.1.fasta file.2.fasta 
```
* Actual code
```
java -Xmx2g -jar /mnt/research/rdp/public/RDPTools/Clustering.jar derep -o derep.gyrB.v1.fa all.seqs.gyrB.v1.ids all.seqs.gyrB.v1.samples fungene_8.8_gyrB_10000_unaligned_protein_seqs.fa fungene_8.8_gyrB_2100_unaligned_protein_seqs.fa fungene_8.8_gyrB_8882_unaligned_protein_seqs.fa 
```
* Results
  * Total sequences: 20982
  * Unique sequences: 4330
* 4,000 is okay, but it still seems a bit high. I also want to check the diversity

## Visualize diversity
* Align sequences using MUSCLE v. 3.8.31
```
muscle -in derep.gyrB.v1.fa -out aligned.derep.gyrB.v1.fa 
```
* Max length: 1617
* Avg length: 820

## Make tree
* Had to submid a job since there are so many sequences (4k)
* Following .qsub was used
```
#!/bin/bash -login
 
### define resources needed:
### walltime - how long you expect the job to run
#PBS -l walltime=24:00:00
 
### nodes:ppn - how many nodes & cores per node (ppn) that you require
#PBS -l nodes=15:ppn=1
 
### mem: amount of memory that the job will need
#PBS -l mem=10gb
 
### you can give your job a name for easier identification
#PBS -N gyrB.raxml
 
### load necessary modules, e.g.
module load raxml/8.0.6
 
### change to the working directory where your code is located
cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/gyrB
 
### call your executable
raxmlHPC -m PROTGAMMAJTTF -p 12345 -s aligned.derep.gyrB.v1.fa -n gyrb.v1
```
