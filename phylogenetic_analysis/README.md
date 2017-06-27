# Phylogenetic analysis of xander-assembled genes
### Taylor Dunivin


## 1. Prepare phylogenies for quality checks
   * Calls on [phylo.sh](https://github.com/ShadeLab/Xander_arsenic/blob/master/phylogenetic_analysis/bin/phylo.sh)
   * To execute, `./phylo.sh GENE CLUST`
   * Should run for each gene with 0.03 as well as 0.01
   * Can be executed from any directory, but working directory should be updated within script per individual    * Inputs: 
      * Gene directory should already exist with `reference_seqs.fa` within it
        * fasta file should include 
          1) seed sequences for gene: from `RDPTools/Xander_assembler/gene_resource/GENE/originaldata/GENE.seeds`
          2) top blast hits for gene: `ncbi.input.txt` -> [entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) -> fasta sequences
          3) a sequence to root by (origin varies) 
      * `coverage.txt`: ALL coverage files from gene of interest; output from xander `search`
      * `final_prot_aligned.fasta` : All aligned protein files for gene of interest; output from xander `search`
   * Outputs: 
      * outputs are flagged with `_short` because they include short (<90% of hmm lenght) sequences
      * `GENE_CLUST_tree_short.nwk`: tree file for iTOL upload  
      * `GENE_CLUST_labels_short.txt`: comma separated labels ultimately for iTOL labeling file
    
## 2. Quality check each gene of interest using phylogenies   
   * [iTOL](http://itol.embl.de/)
      * upload `${GENE}_${CLUST}_tree.nwk`
      * add labels from `${GENE}_${CLUST}_labels.txt` to iTOL labels script
      * [Example file](http://itol.embl.de/help/labels_template.txt)
      * Drag and drop labels onto tree in interface
   * Examine tree for outliers, clustering outside of target sequences
   
## 3. Examine phylogeny of assembled genes
   * Essentially the same as steps 1 and 2 _except_ these scripts are gene specific and remove all sequences that are less than 90% of hmm length
   * Example script [acr3.phylo.sh](https://github.com/ShadeLab/Xander_arsenic/blob/master/phylogenetic_analysis/bin/acr3.phylo.sh)
   * To execute, `./acr3.phylo.sh acr3 CLUST`
   * should repeat for `0.03` and `0.1` once again
   * Inputs are the same as for `phylo.sh`
   * Outputs are the same but do __not__ contain `_short` flag
   * Upload to iTOL and examine diversity ! 
   * Table below shows sequence cutoffs used
   
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
   
## 4. (optional) Group related sequences and examine phylogenies
   * Allows rooting to be completed with other proteins of interest & checks model quality/ specificty
   * Relevant for the following pairs
       * arrA, arxA, aioA
       * arsB, acr3
       * ClassA, ClassB, ClassC
       * AAC3-Ia, AAC6-Ia, AAC6-Ib
       * arsC_glut, arsC_thio
   * Use `group.muscle.sh` for full length sequences and `short.group.muscle.sh` for all sequences
   * To execute, `./group.muscle.sh GROUP GENE1 GENE2 GENE3 CLUST`
       * `GROUP` = name of group; directory for this group should already exist with `reference_seqs.fa` in it
       * `GENE123` = will take up to 3 genes; if you are grouping two genes, simply put `NA` for gene 3
       * `CLUST` = what cluster cutoff you would like to run
       * Aligns using MUSCLE rather than hmmaligner since it takes sequences from up to three HMMs
   * View in iTOL (note: do not need to adjust labels here)
