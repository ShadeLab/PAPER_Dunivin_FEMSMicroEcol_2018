### Taylor Dunivin
## April 5, 2017
---
Short script to merge unaligned protein and nucleotide sequences
Includes dereplication for nucleotide sequences _only_. 
```
#merge sequence files (FunGene can only download 10k at a time)
#do NOT do this for aligned sequences
#this step also makes sure files are named the same (regardless of gene) for automation
cat *unaligned_nucleotide_seqs_v1.fa >merged_unaligned_nucleotide_seqs_v1.fa
cat *unaligned_protein_seqs_v1.fa>framebot.fa

#remove duplicate accession numbers
awk '/^>/{f=!d[$1];d[$1]=1}f' merged_unaligned_nucleotide_seqs_v1.fa >derepaccno.input.fa

#dereplicate nucleotide sequences
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar derep -o nucl.fa derep.all_seqs.ids derep.all_seqs.samples derepaccno.input.fa
```
