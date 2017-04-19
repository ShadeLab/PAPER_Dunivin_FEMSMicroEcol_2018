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
#dereplicate
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar derep -m '#=GC_RF' -o derep.fa all_seqs.ids all_seqs.samples alignment/*.fasta

#calculate distance matrix
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar dmatrix --id-mapping all_seqs.ids --in derep.fa --outfile derep_matrix.bin -l 200  --dist-cutoff 0.1

#cluster
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar cluster --dist-file derep_matrix.bin --id-mapping all_seqs.ids --sample-mapping all_seqs.samples --method complete --outfile all_seq_complete.clust

#export to r format (do this for several clustering distances)
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar cluster_to_Rformat all_seq_complete.clust . 0 0.03

##REPRESENTATIVE SEQUENCES
#merge aligned fasta files
java -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/AlignmentTools/dist/AlignmentTools.jar alignment-merger alignment merged_aligned.fasta
  
#get and rename sequences
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar rep-seqs -c --id-mapping all_seqs.ids --one-rep-per-otu all_seq_complete.clust 0.03 merged_aligned.fasta

#create biom file
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar cluster-to-biom all_seq_complete.clust 0.03 > all_seq_complete.clust.biom

#add info
biom add-metadata -i all_seq_complete.clust.biom -o all_seq_complete.clust.meta.biom --sample-metadata-fp ../meta.txt


java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar derep -f -o all_seq_complete.clust_rep_seqs_modelonly.fasta ids samples all_seq_complete.clust_rep_seqs.fasta
```

Great website for learning phyloseq: https://joey711.github.io/phyloseq/import-data.html

Need to adjust sequence names, but this makes a .nwk tree of final contigs
```
fasttree -nt -gtr < all_seq_complete.clust_rep_seqs_modelonly.fasta > my_expt_tree.nwk
```
