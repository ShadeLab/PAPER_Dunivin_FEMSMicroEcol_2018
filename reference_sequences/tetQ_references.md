# IntI: diverse gene references
### Taylor Dunivin
## March 6, 2017
### Goals: 
* Obtain sets of relatively high quality sequences
* Dereplicate sequences
* Make tree to detemine diversity of sequences
* Assess (and repeat)

## TetQ
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

