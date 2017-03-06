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
