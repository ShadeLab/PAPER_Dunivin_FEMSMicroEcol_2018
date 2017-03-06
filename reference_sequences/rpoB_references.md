# rpoB diverse gene references
### Taylor Dunivin
## March 6, 2017
### Goals: 
* Obtain a set of relatively high quality RpoB sequences
* Dereplicate sequences
* Make tree to detemine diversity of sequences
* Assess (and repeat)

## Downloading from FunGene
* RpoB FunGene database already exists
  * authored by Scott Santos/Howard Ochman
  * 146173 existing sequences
* Need to narrow down based on
  1. HMM coverage - typically don't go below 80%
  2. Length (aa) - RpoB is 1193 aa
  3. Score - based on where there is a large score drop off

| Protein | Test # | min aa | min HMM coverage (%) | min HMM score |
| --------- | ----- | ---------- | --------- | :-----: |
| RpoB | 1 | 1175 | 30 | 1000 |

* There is a low HMM coverage here because some of the seed sequences are >4000aa sequences, so even full length options (1193aa) have low HMM coverages
* This leaves 93482 sequences for derep
* Derep sequences
```
java -Xmx2g -jar /mnt/research/rdp/public/RDPTools/Clustering.jar derep -o derep.rpoB.v1.fa all.seqs.rpoB.v1.ids all.seqs.rpoB.v1.samples fungene_8.8_rpoB_10000_unaligned_protein_seqs11-15.fa  fungene_8.8_rpoB_10000_unaligned_protein_seqs36-40.fa fungene_8.8_rpoB_10000_unaligned_protein_seqs16-20.fa  fungene_8.8_rpoB_10000_unaligned_protein_seqs41-45.fa fungene_8.8_rpoB_10000_unaligned_protein_seqs21-25.fa  fungene_8.8_rpoB_10000_unaligned_protein_seqs6-10.fa fungene_8.8_rpoB_10000_unaligned_protein_seqs26-30.fa  fungene_8.8_rpoB_3482_unaligned_protein_seqs46-47.fa fungene_8.8_rpoB_10000_unaligned_protein_seqs31-35.fa  fungene_8.8_rpoB_9025_unaligned_protein_seqs1-5.fa
```

