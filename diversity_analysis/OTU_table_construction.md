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
```

The RDP has a script set up to specifically make OTU tables for R from Xander outputs. Use this in the following steps to get OTUs and representative sequences
```
##CLUSTERING
../get_OTUabundance.sh final_coverage.txt /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/OTUabundances/GENE 0 0.1 alignment/*

#outputs:
#ids - gives OTU number and contig name
#samples - gives contig name and site of origin
#rformat* - gives R formatted OTU table for each specified distance cutoff
#derep.fa - aligned fasta file of all dereplicated sequences with contig name as heading
#complete.clust - gives cluster number with corresponding site, cluster count, and contig name

##REPRESENTATIVE SEQUENCES
#get and rename sequences
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar rep-seqs -c --id-mapping ids --one-rep-per-otu complete.clust 0.03 derep.fa

#outputs:
#complete.clust_rep_seqs.fasta - gives aligned file with cluster (OTU) number and represenative sequences at specified distance cutoff
```

The representative sequences unfortunately have a nomenclature issue where they are called "cluster" while the others are called "OTU." We can simple replace cluster with otu
```
sed -i 's/cluster_/OTU_/g' complete.clust_rep_seqs.fasta
```

Now that sequence names have been adjusted, but we can make a .nwk tree of representative sequences of final OTUs
```
module load GNU/4.4.5
module load FastTree/2.1.7
FastTree complete.clust_rep_seqs_modelonly.fasta > gene_dist_tree.nwk
```



Great website for learning phyloseq: https://joey711.github.io/phyloseq/import-data.html
