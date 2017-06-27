#!/bin/bash

#requirements: gene group name, gene1 (will be used for model), gene 2, gene 3 (optional), cluster cutoff

#you have to type the gene
GROUP=$1
GENE1=$2
GENE2=$3
GENE3=$4
CLUST=$5


#change to working directory
cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances

#make directory for group
mkdir ${GROUP}
cd ${GROUP}

#get all coverage files 
#copy all files to database
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/${GENE1}/complete.clust_rep_seqs_${CLUST}_unaligned_short.fasta ${GENE1}_complete.clust_full_rep_seqs_${CLUST}_unaligned_short.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/${GENE2}/complete.clust_rep_seqs_${CLUST}_unaligned_short.fasta ${GENE2}_complete.clust_full_rep_seqs_${CLUST}_unaligned_short.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/${GENE3}/complete.clust_rep_seqs_${CLUST}_unaligned_short.fasta ${GENE3}_complete.clust_full_rep_seqs_${CLUST}_unaligned_short.fasta

#change OTU to gene name
sed -i 's/OTU_/${GENE1}_/g' ${GENE1}_complete.clust_full_rep_seqs_${CLUST}_unaligned_short.fasta
sed -i 's/OTU_/${GENE2}_/g' ${GENE2}_complete.clust_full_rep_seqs_${CLUST}_unaligned_short.fasta
sed -i 's/OTU_/${GENE3}_/g' ${GENE3}_complete.clust_full_rep_seqs_${CLUST}_unaligned_short.fasta

#add seed data to sequence file
#Note: first would need to uploade sequences to HPCC
cat *complete.clust_full_rep_seqs_${CLUST}_unaligned_short.fasta reference_seqs.fa > complete.clust_full_rep_seqs_${CLUST}_unaligned_seeds_short.fasta

#align sequences
module load muscle
muscle -in complete.clust_full_rep_seqs_${CLUST}_unaligned_seeds_short.fasta -out ${GROUP}_alignment_muscle_${CLUST}_short.fa

#make tree
module load FastTree
FastTree ${GROUP}_alignment_muscle_${CLUST}_short.fa > ${GROUP}_muscle_${CLUST}_short.nwk


#extract labels from fasta file
grep "^>" ${GROUP}_alignment_muscle_${CLUST}_short.fa > ${GROUP}_alignment_muscle_${CLUST}_labels_short.txt

#edit file name for iTOL
sed 's/^........../,/' ${GROUP}_alignment_muscle_${CLUST}_labels_short.txt > ${GROUP}_${CLUST}_labels_short.txt
sed 's/, /,/' ${GROUP}_${CLUST}_labels_short.txt > ${GROUP}_${CLUST}_labels_short.n.txt

#add numbers to each line (seq_name)
nl -w 1 -s "" ${GROUP}_${CLUST}_labels_short.n.txt > ${GROUP}_${CLUST}_labels_short.txt

#delete unnecessary files
rm ${GROUP}_${CLUST}_labels_short.n.txt

