# TetQ: diverse gene references
### Taylor Dunivin
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
