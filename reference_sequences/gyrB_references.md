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

```
java -Xmx2g -jar /path/to/Clustering.jar derep -o derep.gyrB.fa all.seqs.gyrB.ids all.seqs.gyrB.samples file.1.fasta file.2.fasta 
```
