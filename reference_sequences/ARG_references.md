# Antibiotic resistance genes: diverse gene references
### Taylor Dunivin
## March 6, 2017
### Goals: 
* Obtain sets of relatively high quality sequences
* Dereplicate sequences
* Make tree to detemine diversity of sequences
* Assess (and repeat)

## IntI
### Downloading from FunGene
* IntI FunGene database already exists
  * authored by Carlos Rodriguez-Minguela
  * 137607 existing sequences
* Need to narrow down based on
  1. HMM coverage - typically don't go below 80%
  2. Length (aa) - IntI is typically 337 aa
  3. Score - based on where there is a large score drop off

| Protein | Test # | min aa | min HMM coverage (%) | min HMM score |
| --------- | ----- | ---------- | --------- | :-----: |
| IntI | 1 | 315 | 80 | 90 |

* This resulted in 5444 sequences, which is too many.
* I will see how dereplication impacts this
```
java -Xmx2g -jar /mnt/research/rdp/public/RDPTools/Clustering.jar derep -o derep.intI.v1.fa all.seqs.intI.v1.ids all.seqs.intI.v1.samples fungene_8.8_intI_5444_unaligned_protein_seqs.fa
```
* Total unique sequences: 1169
* I also need to remove distant sequences (some might not be intI)

### Aligning sequences

## TetQ
### Gather sequences
* TetQ FunGene database already exists
 * authored by
 * ## sequences exist
 
| Protein | Test # | min aa | min HMM coverage (%) | min HMM score |
| --------- | ----- | ---------- | --------- | :-----: |
| TetQ | 1 | 620 | 80 | 900 |

* Results in 233 sequences

### Dereplicate sequences
```
java -Xmx2g -jar /mnt/research/rdp/public/RDPTools/Clustering.jar derep -o derep.tetQ.v1.fa all.seqs.tetQ.v1.ids all.seqs.tetQ.v1.samples fungene_8.8_tetQ_233_unaligned_protein_seqs.fa
```
=56 sequences

### Align sequences
```
muscle -in derep.tetQ.v1.fa -out aligned.derep.tetQ.v1.fa
```
Longest: 670
Avg: 641
