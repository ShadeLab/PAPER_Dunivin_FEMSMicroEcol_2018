#!/bin/bash

#you have to type the gene
GENE=$1
CLUST=$2

#change to working directory
cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/${GENE}

#make alignment directory
mkdir alignment
cd alignment

#copy all aligned protein files 
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen01/k45/${GENE}/cluster/cen01_${GENE}_45_final_prot_aligned.fasta Cen01.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen03/k45/${GENE}/cluster/cen03_${GENE}_45_final_prot_aligned.fasta Cen03.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen04/k45/${GENE}/cluster/cen04_${GENE}_45_final_prot_aligned.fasta Cen04.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen05/k45/${GENE}/cluster/cen05_${GENE}_45_final_prot_aligned.fasta Cen05.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen06/k45/${GENE}/cluster/cen06_${GENE}_45_final_prot_aligned.fasta Cen06.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen07/k45/${GENE}/cluster/cen07_${GENE}_45_final_prot_aligned.fasta Cen07.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen10/k45/${GENE}/cluster/cen10_${GENE}_45_final_prot_aligned.fasta Cen10.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen12/k45/${GENE}/cluster/cen12_${GENE}_45_final_prot_aligned.fasta Cen12.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen13/k45/${GENE}/cluster/cen13_${GENE}_45_final_prot_aligned.fasta Cen13.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen14/k45/${GENE}/cluster/cen14_${GENE}_45_final_prot_aligned.fasta Cen14.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen15/k45/${GENE}/cluster/cen15_${GENE}_45_final_prot_aligned.fasta Cen15.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen16/k45/${GENE}/cluster/cen16_${GENE}_45_final_prot_aligned.fasta Cen16.fasta
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen17/k45/${GENE}/cluster/cen17_${GENE}_45_final_prot_aligned.fasta Cen17.fasta

#move up a directory
cd ..

#get all coverage files 
#copy all files to database
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen01/k45/${GENE}/cluster/cen01_${GENE}_45_coverage.txt cen01_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen03/k45/${GENE}/cluster/cen03_${GENE}_45_coverage.txt cen03_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen04/k45/${GENE}/cluster/cen04_${GENE}_45_coverage.txt cen04_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen05/k45/${GENE}/cluster/cen05_${GENE}_45_coverage.txt cen05_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen06/k45/${GENE}/cluster/cen06_${GENE}_45_coverage.txt cen06_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen07/k45/${GENE}/cluster/cen07_${GENE}_45_coverage.txt cen07_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen10/k45/${GENE}/cluster/cen10_${GENE}_45_coverage.txt cen10_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen12/k45/${GENE}/cluster/cen12_${GENE}_45_coverage.txt cen12_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen13/k45/${GENE}/cluster/cen13_${GENE}_45_coverage.txt cen13_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen14/k45/${GENE}/cluster/cen14_${GENE}_45_coverage.txt cen14_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen15/k45/${GENE}/cluster/cen15_${GENE}_45_coverage.txt cen15_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen16/k45/${GENE}/cluster/cen16_${GENE}_45_coverage.txt cen16_${GENE}_45_coverage.txt
cp /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/cen17/k45/${GENE}/cluster/cen17_${GENE}_45_coverage.txt cen17_${GENE}_45_coverage.txt

#merge coverage.txt files
cat *45_coverage.txt > final_coverage.txt

##CLUSTERING
../bin/./get_OTUabundance.sh final_coverage.txt /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/${GENE} 0 ${CLUST} alignment/*

#rename rformat of interest
mv rformat_dist_${CLUST}.txt ${GENE}_rformat_dist_${CLUST}.txt

#get and rename sequences
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar rep-seqs -c --id-mapping ids --one-rep-per-otu complete.clust ${CLUST} derep.fa

#rename complete.clust_rep_seqs.fasta to match distance
mv complete.clust_rep_seqs.fasta complete.clust_rep_seqs_${CLUST}.fasta 

#rename files based on OTU, not cluster for specified distance
sed -i 's/[0-9]\{1,\}/0000000&/g;s/0*\([0-9]\{4,\}\)/\1/g' complete.clust_rep_seqs_${CLUST}.fasta
sed -i 's/cluster_/OTU_/g' complete.clust_rep_seqs_${CLUST}.fasta 

#unalign sequences
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar to-unaligned-fasta complete.clust_rep_seqs_${CLUST}.fasta > complete.clust_rep_seqs_${CLUST}_unaligned.fasta

#remove sequences below desired length
#####NUMBER SHOULD BE 0.9*HMM LENGTH AA#######
#####CHANGE FOR EACH GENE OF INTEREST!!#######
../bin/bioawk/./bioawk -c fastx '{ if(length($seq) > 327) { print ">"$name; print $seq }}' complete.clust_rep_seqs_${CLUST}_unaligned.fasta > complete.clust_full_rep_seqs_${CLUST}_unaligned.fasta

#add seed data to sequence file
#Note: first would need to uploade sequences to HPCC
cat complete.clust_full_rep_seqs_${CLUST}_unaligned.fasta reference_seqs.fa > complete.clust_full_rep_seqs_${CLUST}_unaligned_seeds.fasta

#make tree for visualization:
#align file with all sequences
module load HMMER/3.1b2
hmmalign --amino --outformat SELEX -o ${GENE}_alignment_${CLUST}.selex ../../analysis/RDPTools/Xander_assembler/gene_resource/${GENE}/originaldata/${GENE}.hmm complete.clust_full_rep_seqs_${CLUST}_unaligned_seeds.fasta

#convert alignment from selex format to aligned fasta (xmfa)
module load Bioperl/1.6.923
../bin/./convertAlignment.pl -i ${GENE}_alignment_${CLUST}.selex -o ${GENE}_alignment_${CLUST}.fa -f xmfa -g selex

#remove last line of aligned seqs (= sign)
dd if=/dev/null of=${GENE}_alignment_${CLUST}.fa bs=1 seek=$(echo $(stat --format=%s ${GENE}_alignment_${CLUST}.fa ) - $( tail -n1 ${GENE}_alignment_${CLUST}.fa | wc -c) | bc )

#make two separate trees
module load GNU/4.4.5
module load FastTree/2.1.7
FastTree ${GENE}_alignment_${CLUST}.fa > ${GENE}_${CLUST}_tree.nwk

#extract labels from fasta file
grep "^>" ${GENE}_alignment_${CLUST}.fa > ${GENE}_alignment_${CLUST}_labels.txt

#edit file name for iTOL
sed 's/^........../,/' ${GENE}_alignment_${CLUST}_labels.txt > ${GENE}_${CLUST}_labels.txt
sed 's/, /,/' ${GENE}_${CLUST}_labels.txt > ${GENE}_${CLUST}_labels.n.txt

#add numbers to each line (seq_name)
nl -w 1 -s "" ${GENE}_${CLUST}_labels.n.txt > ${GENE}_${CLUST}_labels.txt

#delete unnecessary files
rm rformat_dist_*
rm ${GENE}_${CLUST}_labels.n.txt
rm *_${GENE}_45_coverage.txt
