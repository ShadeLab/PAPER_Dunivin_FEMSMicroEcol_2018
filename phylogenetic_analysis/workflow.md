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
