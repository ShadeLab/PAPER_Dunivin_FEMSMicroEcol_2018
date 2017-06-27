#!/bin/bash

#you have to type the gene
GENE=$1

#change to working directory
cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/blast/${GENE}

#copy all files to database
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen01/k45/${GENE}/cluster/cen01_${GENE}_45_final_prot.fasta cen01_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen03/k45/${GENE}/cluster/cen03_${GENE}_45_final_prot.fasta cen03_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen04/k45/${GENE}/cluster/cen04_${GENE}_45_final_prot.fasta cen04_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen05/k45/${GENE}/cluster/cen05_${GENE}_45_final_prot.fasta cen05_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen06/k45/${GENE}/cluster/cen06_${GENE}_45_final_prot.fasta cen06_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen07/k45/${GENE}/cluster/cen07_${GENE}_45_final_prot.fasta cen07_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen10/k45/${GENE}/cluster/cen10_${GENE}_45_final_prot.fasta cen10_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen12/k45/${GENE}/cluster/cen12_${GENE}_45_final_prot.fasta cen12_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen13/k45/${GENE}/cluster/cen13_${GENE}_45_final_prot.fasta cen13_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen14/k45/${GENE}/cluster/cen14_${GENE}_45_final_prot.fasta cen14_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen15/k45/${GENE}/cluster/cen15_${GENE}_45_final_prot.fasta cen15_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen16/k45/${GENE}/cluster/cen16_${GENE}_45_final_prot.fasta cen16_${GENE}_45_final_prot.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen17/k45/${GENE}/cluster/cen17_${GENE}_45_final_prot.fasta cen17_${GENE}_45_final_prot.fasta

#blast nr database!
module load BLAST/2.2.26 
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB
cat *_final_prot.fasta >final_prot.fasta
blastall -d nr -i final_prot.fasta -p blastp -o blast.txt -b 1 -v 1 -e 1e-6 -a 8
grep '^>' blast.txt > descriptor.blast.txt

#get gene descriptions
cat descriptor.blast.txt | awk -F '[|]' '{print $3}' > gene.descriptor.txt

#count occurrences of gene descriptions
sort gene.descriptor.txt | uniq -c > gene.descriptor.final.txt

#remove unnecessary files
rm gene.descriptor.txt

#get and count occurrences of accession numbers
cat descriptor.blast.txt | awk -F '[|]' '{print $2}' > accno.txt

#count occurrences of gene descriptions
sort accno.txt | uniq -c > accno.final.txt

cat accno.final.txt | awk '{print $2}' > ncbi.input.txt

#remove unnecessary files
rm accno.txt
rm *_${GENE}_45_final_prot.fasta
