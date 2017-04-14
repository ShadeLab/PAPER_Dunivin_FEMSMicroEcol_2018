```
##DATA PREP
#move coverage files of the same gene from all centralia sites to a separate analysis folder
cp *coverage.txt /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/rplB

#move aligned protein files of the same gene from all centralia sites to a separate analysis folder
cp *_final_prot_aligned.fasta /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/rplB/alignment

#move to gene directory
cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/GENE

#merge coverage.txt files
cat *coverage.txt >final_coverage.txt

##CLUSTERING
#get OTUabundance
#input: all coverage information output directory min_distance max_distance aligned_prot.fasta
./../get_OTUabundance.sh final_coverage.txt /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/rplB 0 0.05 alignment/*

#output: ids, samples, derep.fa, complete.clust, rformat_dist_0.0*

##REPRESENTATIVE SEQUENCES
#merge aligned fasta files
java -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/AlignmentTools/dist/AlignmentTools.jar alignment-merger alignment merged_aligned.fasta

#get and rename sequences
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar rep-seqs -c --id-mapping ids --one-rep-per-otu complete.clust 0.03 merged_aligned.fasta

#create biom file
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar cluster-to-biom complete.clust 0.03 > all_seq_complete.clust.biom

#add info
biom add-metadata -i all_seq_complete.clust.biom -o all_seq_complete.clust.meta.biom --sample-metadata-fp ../meta.txt

biom add-metadata -i all_seq_complete.clust.biom -o all_seq_complete.clust.meta.biom --sample-metadata-fp ../meta.txt

java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar derep -f -o all_seq_complete.clust_rep_seqs_modelonly.fasta ids samples all_seq_complete.clust_rep_seqs.fasta

https://joey711.github.io/phyloseq/import-data.html

fasttree -nt -gtr < all_seq_complete.clust_rep_seqs_modelonly.fasta > my_expt_tree.nwk

http://www.jennajacobs.org/R/rarefaction.html
```

```
#might need?
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar derep -m '#=GC_RF' -o derep.fa all_seqs.ids all_seqs.samples alignment/*.fasta

java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar dmatrix --id-mapping all_seqs.ids --in derep.fa --outfile derep_matrix.bin   -l 200 --dist-cutoff 0.1

java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar cluster --dist-file derep_matrix.bin --id-mapping all_seqs.ids --sample-mapping all_seqs.samples --method complete --outfile all_seq_complete.clust

java -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/AlignmentTools/dist/AlignmentTools.jar alignment-merger alignment merged_aligned.fasta

java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar rep-seqs -c --id-mapping all_seqs.ids --one-rep-per-otu   all_seq_complete.clust 0.03 merged_aligned.fasta



java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/classifier.jar classify -c 0.5 -f biom -m all_seq_complete.clust.biom -d sam.data.txt -o  all_seq_complete.clust_classified.biom complete.clust_rep_seqs.fasta
```
