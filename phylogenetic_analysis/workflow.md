# Workflow/ notes for phylogenetic analysis
### Taylor Dunivn

# Table of contents
* [June 9, 2017]()
* [June 20, 2017]()


# June 9, 2017
```
#add together all coverage information
cat *coverage.txt >final_coverage.txt

#cluster with abundance information
../get_OTUabundance.sh final_coverage.txt /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/aioA 0 0.5 alignment/*


##get representative sequences for BOTH 0.1 and 0.2 sequence distance
#0.1
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar rep-seqs -c --id-mapping ids --one-rep-per-otu complete.clust 0.1 derep.fa
#rename complete.clust_rep_seqs.fasta to match distance
mv complete.clust_rep_seqs.fasta complete.clust_rep_seqs_0.1.fasta 

#0.2
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar rep-seqs -c --id-mapping ids --one-rep-per-otu complete.clust 0.2 derep.fa
#rename complete.clust_rep_seqs.fasta to match distance
mv complete.clust_rep_seqs.fasta complete.clust_rep_seqs_0.2.fasta 

#rename files based on OTU, not cluster for both 0.1 and 0.2 distances
sed -i 's/[0-9]\{1,\}/0000000&/g;s/0*\([0-9]\{2,\}\)/\1/g' complete.clust_rep_seqs_0.1.fasta 
sed -i 's/cluster_/OTU_/g' complete.clust_rep_seqs_0.1.fasta 

sed -i 's/[0-9]\{1,\}/0000000&/g;s/0*\([0-9]\{2,\}\)/\1/g' complete.clust_rep_seqs_0.2.fasta
sed -i 's/cluster_/OTU_/g' complete.clust_rep_seqs_0.2.fasta

#unalign sequences
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar to-unaligned-fasta complete.clust_rep_seqs_0.1.fasta > complete.clust_rep_seqs_0.1_unaligned.fasta
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar to-unaligned-fasta complete.clust_rep_seqs_0.2.fasta > complete.clust_rep_seqs_0.2_unaligned.fasta

#add seed data to sequence file
#Note: first would need to uploade sequences to HPCC
cat complete.clust_rep_seqs_0.1_unaligned.fasta aioA_reference_seqs.fa > complete.clust_rep_seqs_0.1_unaligned_seeds.fasta
cat complete.clust_rep_seqs_0.2_unaligned.fasta aioA_reference_seqs.fa > complete.clust_rep_seqs_0.2_unaligned_seeds.fasta

#align file with all sequences
module load HMMER/3.1b2
hmmalign --amino --outformat SELEX -o aioA_alignment_0.1.selex ../../analysis/RDPTools/Xander_assembler/gene_resource/aioA/originaldata/aioA.hmm complete.clust_rep_seqs_0.1_unaligned_seeds.fasta
hmmalign --amino --outformat SELEX -o aioA_alignment_0.2.selex ../../analysis/RDPTools/Xander_assembler/gene_resource/aioA/originaldata/aioA.hmm complete.clust_rep_seqs_0.2_unaligned_seeds.fasta

#convert alignment from selex format to aligned fasta (xmfa)
module load Bioperl/1.6.923
.././convertAlignment.pl -i aioA_alignment_0.1.selex -o aioA_alignment_0.1.fa -f xmfa -g selex
#need to remove = at bottom
.././convertAlignment.pl -i aioA_alignment_0.2.selex -o aioA_alignment_0.2.fa -f xmfa -g selex
#need to remove = at bottom

#make two separate trees
module load GNU/4.4.5
module load FastTree/2.1.7
FastTree aioA_alignment_0.1.fa > aioA_0.1_tree.nwk
FastTree aioA_alignment_0.2.fa > aioA_0.2_tree.nwk

#extract labels from fasta file
grep "^>" aioA_alignment_0.1.fa > aioA_alignment_0.1_labels.txt
grep "^>" aioA_alignment_0.2.fa > aioA_alignment_0.2_labels.txt

#edit file name for iTOL
sed 's/^........../,/' aioA_alignment_0.1_labels.txt > aioA_0.1_labels.txt
sed 's/, /,/' aioA_0.1_labels.txt > aioA_0.1_labels.n.txt

sed 's/^........../,/' aioA_alignment_0.2_labels.txt > aioA_0.2_labels.txt
sed 's/, /,/' aioA_0.2_labels.txt > aioA_0.2_labels.n.txt

#add numbers to each line (seq_name)
nl -w 1 -s "" aioA_0.1_labels.n.txt >aioA_0.1_labels.txt
nl -w 1 -s "" aioA_0.2_labels.n.txt >aioA_0.2_labels.txt

#export to computer
#aioA_0.1_labels.txt
#aioA_0.2_labels.txt
#aioA_0.1_tree.nwk
```

# June 20, 2017
* Goals
  * Write scripts to cluster sequences at specified cutoffs, unalign, incorporate reference sequences, make maximum likelihood tree, and extract OTU labels (for iTOL format)
  * Add script version to cut off size of sequences included to >= 90% hmm length (more appropriate to do phylogenetic analysis on full length sequences: Edwardson, C. F., & Hollibaugh, J. T. (2017). Metatranscriptomic analysis of prokaryotic communities active in sulfur and arsenic cycling in Mono Lake, California, USA. The ISME Journal. https://doi.org/10.1038/ismej.2017.80)
  * Add script to group relevant sequences for properly booted tree/ comparison

* Script to perform initial analysis
  * Stored in ```/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/bin```
  * Titled: ```phylo.sh```
  * Pre-script activities
    * create directory of gene name
    * add file `reference_seqs.fa` 
      * Includes seed sequences for gene(s) at play
      * Includes top blast hits (see above; output `ncbi.input.txt` from `blast.summary.pl`) 
  * Relevant outputs: 
    * `${GENE}_${CLUST}_labels_short.txt`: Labels of all sequences (short) for incorporating into iTOL trees
    * `${GENE}_${CLUST}_tree_short.nwk`: Maximum likelihood tree of all sequences (short) for iTOL tree
  * Script
```
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
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar to-unaligned-fasta complete.clust_rep_seqs_${CLUST}.fasta > complete.clust_rep_seqs_${CLUST}_unaligned_short.fasta

#add seed data to sequence file
#Note: first would need to uploade sequences to HPCC
cat complete.clust_full_rep_seqs_${CLUST}_unaligned_short.fasta reference_seqs.fa > complete.clust_full_rep_seqs_${CLUST}_unaligned_seeds_short.fasta

#make tree for visualization:
#align file with all sequences
module load HMMER/3.1b2
hmmalign --amino --outformat SELEX -o ${GENE}_alignment_${CLUST}_short.selex ../../analysis/RDPTools/Xander_assembler/gene_resource/${GENE}/originaldata/${GENE}.hmm complete.clust_full_rep_seqs_${CLUST}_unaligned_seeds_short.fasta

#convert alignment from selex format to aligned fasta (xmfa)
module load Bioperl/1.6.923
../bin/./convertAlignment.pl -i ${GENE}_alignment_${CLUST}_short.selex -o ${GENE}_alignment_${CLUST}_short.fa -f xmfa -g selex

#remove last line of aligned seqs (= sign)
dd if=/dev/null of=${GENE}_alignment_${CLUST}_short.fa bs=1 seek=$(echo $(stat --format=%s ${GENE}_alignment_${CLUST}_short.fa ) - $( tail -n1 ${GENE}_alignment_${CLUST}_short.fa | wc -c) | bc )

#make two separate trees
module load GNU/4.4.5
module load FastTree/2.1.7
FastTree ${GENE}_alignment_${CLUST}_short.fa > ${GENE}_${CLUST}_tree_short.nwk

#extract labels from fasta file
grep "^>" ${GENE}_alignment_${CLUST}_short.fa > ${GENE}_alignment_${CLUST}_labels_short.txt

#edit file name for iTOL
sed 's/^........../,/' ${GENE}_alignment_${CLUST}_labels_short.txt > ${GENE}_${CLUST}_labels_short.txt
sed 's/, /,/' ${GENE}_${CLUST}_labels_short.txt > ${GENE}_${CLUST}_labels_short.n.txt

#add numbers to each line (seq_name)
nl -w 1 -s "" ${GENE}_${CLUST}_labels_short.n.txt > ${GENE}_${CLUST}_labels_short.txt

#delete unnecessary files
rm rformat_dist_*
rm ${GENE}_${CLUST}_labels_short.n.txt
rm *_${GENE}_45_coverage.txt
```
 
* Script to perform __gene-length-specific__ analysis
  * Stored in ```/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/bin```
  * Titled: ```<gene>.phylo.sh```
  * Pre-script activities
    * create directory of gene name
    * add file `reference_seqs.fa` 
      * Includes seed sequences for gene(s) at play
      * Includes top blast hits (see above; output `ncbi.input.txt` from `blast.summary.pl`) 
  * Relevant outputs: 
    * `${GENE}_${CLUST}_labels.txt`: Labels of all sequences (short) for incorporating into iTOL trees
    * `${GENE}_${CLUST}_tree.nwk`: Maximum likelihood tree of all sequences (short) for iTOL tree
  * Script (acr3 example)
  
```
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
```

| Gene | HMM length | 90% length | Comments |
| ---- | ---------- | ---------- | -------- |
| acr3 | 363 | 326  |  | 
| aioA | 837 | 753 |  | 
| arsB | 430 | 387 |  | 
| arsCglut | 117 | 105 |  | 
| arsCthio | 133 | 119 |  | 
| arsD | 123 | 110 |  | 
| arsM | 269 | 242 |  | 
| arxA | 831 | 747 |  | 
| arrA | 846 | 761 |  | 
| tolC | 439 | 395 |  | 
| vanZ | 161 | 144 |  | 
| vanX | 202 | 181 |  | 
| vanH | 323 | 290 |  | 
| vanA | 343 | 308 | No sequences >= 90% | 
| tetX | 398 | 358 |  | 
| tetW | 639 | 575 |  | 
| tetA | 406 | 365 |  | 
| sul2 | 272 | 244 | No sequences >= 90% | 
| rplB | 277 | 249 |  | 
| intI | 319 | 287 | No sequences >= 90% | 
| CEP | 416 | 374 |  | 
| CAT | 217 | 195 |  | 
| arsA | 346 | 311 |  | 
| adeB | 1035 | 931 |  | 
| AAC6-Ib | 184 | 165 | No sequences >= 90% | 
| AAC6-Ia | 185 | 166 | No sequences >= 90% | 
| AAC3-Ia | 154 | 138 | No sequences >= 90% | 
| ClassA | 294 | 264 |  | 
| ClassB | 264 | 237 |  | 
| ClassC | 385 | 346 | No sequences >= 90% | 
  
  
  
  
  
  
  
  
  
