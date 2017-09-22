# Phylogenetic analysis of xander-assembled genes
### Taylor Dunivin
For analysis in iTOL
   * Need to extract labels and convert to iTOL format -> [names_iTOL.sh](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/blob/master/assembly_assessments/bin/names_iTOL.sh)
      * Inputs
        * `complete.clust_full_rep_seqs_${CLUST}_unaligned_seeds.fasta`: from phylo.hmm.align.sh script
        * `iTOL.txt`: iTOL header
    * To include relative abundance information, use this R script (or use code as R funciton) [tree_setup.R](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/blob/master/phylogenetic_analysis/tree_setup.R) and adjust colors, etc in iTOL forms. 
     
