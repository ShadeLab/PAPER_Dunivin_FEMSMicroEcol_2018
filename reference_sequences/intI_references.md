# IntI: diverse gene references
### Taylor Dunivin
## March 6, 2017
### Goals: 
* Obtain sets of relatively high quality sequences
* Dereplicate sequences
* Make tree to detemine diversity of sequences
* Assess (and repeat)

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
Since there are not many sequences, I did not need to submit a job.
```
module load MUSCLE/3.8.31
muscle -in derep.intI.v1.fa -out aligned.derep.intI.v1.fa
```
* Longest sequence (aa): 594
* Avg length (aa): 346

## March 13, 2017
### Goals:
* Make trees
* Edit diverse gene ref lists

### Make tree
Since there are not many sequences, I did not need to submit a job.
```
module load raxml/8.0.6
raxmlHPC -m PROTGAMMAJTTF -p 12345 -s aligned.derep.intI.v1.fa -n intI.v1
```

I visualized the trees using iTol
* ![Parsimony tree](https://github.com/ShadeLab/Xander_arsenic/blob/master/images/RAxML.parsimony.intI.png)
 * The tree will not show online since it is too large
 * Git desktop will display it, and downloading it will also show it
 * The tree looks good with no sequence removals apparently necessary
 
 ## March 15, 2017
 ### Goals:
 * Re-download intI sequences with GI instead of accno
 * Also download aligned versions (rather than aligning sequences once they're derep)
 * Derep unaligned protein sequences
 * Derep aligned protein sequences AND unaligned nucleotide sequences (based on derep unaligned prot)
 
 ### Download new sequences
 * Used same filter parameters as before (See above)
 * ```fungene_8.8_intI_9420_unaligned_protein_seqs_v2.fa```
 * ```fungene_8.8_intI_9418_unaligned_nucleotide_seqs_v2.fa```
 * ```fungene_8.8_intI_9421_aligned_protein_seqs_v2.fa```
 
 ### Derep unaligned protein sequences
 ```
 java -Xmx2g -jar /mnt/research/rdp/public/RDPTools/Clustering.jar derep -o derep.intI.v2.fa all.seqs.intI.v2.ids all.seqs.intI.v2.samples fungene_8.8_intI_9420_unaligned_protein_seqs_v2.fa
 ```
 
 * Due to duplicates (because of using GI numbers) one replicate of the following sequences was deleted
  * 151363173
  * 120322793
  * 312940827
  * 1044515116
  * 431827765
  * 239918938
  * 332330813
  * 508120469
  * 310764731
  * 90823168
  * 661563325
  * 119765642
  * 662231943
  * 52214156
  * 1045271457
  * 939199104
  * 161158851
  * 149947715
  * 520999024
  * 945279640
  * 293612874
  * 343923153
  * 946699795
  * 1059051276
  * 223587976
 ...
 * Actually this was not completed. This is not going to be an efficient way to analyze these data. 
 * Instead I have downloaded ftp://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz which contains gene accession numbers in relation to GI numbers
 * I will need to figure out a way to derep with accno, translate accno to GI, translate GI to nucleotide accno, then derep nucleotide based on derep protein data.
  
 
