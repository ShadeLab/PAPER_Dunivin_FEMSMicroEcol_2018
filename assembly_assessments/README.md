# Assembly assessments
### Taylor Dunivin


This directory contains scripts and notes for assessment of Xander's gene targeted assembly. 

## 1. Make blast databases for each gene of interest
   * Need blast databases from reference sequences (used in xander `find` and `search`)
   * Run [makeblastdb.sh](https://github.com/ShadeLab/Xander_arsenic/blob/master/assembly_assessments/bin/makeblastdb.sh)
        * `./makeblastdb.sh GENE`
        * Script can be run in any directory, but working directiry should be changed within the script for individual users
        * Input: `RDPTools/Xander_assembler/gene_resource/${GENE}/originaldata/nucl.fa`
        * Output: gene directory
          * `GENE_database.nhr`
          * `GENE_database.nin`
          * `GENE_database.nsq`

## 2. Assembly assessments
   * For every gene in each sample, assess the following
        * Examine kmer abundance distribution
        * Calculate length (nucl) statistics, % identity statistics, and # at 99%
        * Calculate e-values for all sequences (make sure >10^-5)
   * Run [assembly_assessments.sh](https://github.com/ShadeLab/Xander_arsenic/blob/master/assembly_assessments/bin/xander_assessments.sh)
        * Script mostly authored by ACRES REU student Susannah Yeh
        * `./assembly_assessments.sh GENE SAMPLE`
        * Script can be run in any directory, but working directiry should be changed within the script for individual users
        * Calls on `assembly_assessmentR.R` and `get_gc_counts.pl`
       * Requirements:
          * ```GNU/4.4.5```, ```GNU/4.9```
          * ```Gblastn/2.28```
          * ```OpenMPI/1.10.0```
          * ```R/3.3.0```
          * ```perl/5.24.1```
       * Inputs: 
          * ```final_nucl.fasta```: Final nucleotide sequences from xander `search` output
          * ```match_reads.fa```: Matched reads from xander `search` output
          * ```framebot.txt```: Alignment statistics & alignments from xander `search` output
          * ```/blastdatabases/<gene>_database```: Blast database of gene of interest output from ```makeblastdb.sh```
       * Outputs:
          * GENE/
            * ```stats.txt```: # of protein clusters (99% identity); average, median, min, max length aa; max, min, average % identity
            * ```kmerabundancedist.png```: plot of kmer abundance distribution
            * ```e.values.txt ```: # of sequences with e values < 1e-5; min, max, average, and median e-values
            * ```readssummary.txt```: # unique reads, # total reads
          * GC/
            * ```gc_out.txt```: lists %GC, total bp, #G, #C, #A, #T per assembled contig
            
## 3. Assessments by gene
   * This step accomplishes the following goals
     * Be certain that when considering nr database, the gene of interest is still the top hit
     * Add to confidence in assembly assessments
     * Determine whether/ how to troubleshoot low quality results from specific genes
   * Calls on [blast.summary.pl](https://github.com/ShadeLab/Xander_arsenic/blob/master/assembly_assessments/bin/blast.summary.pl)
      * Script was written in part by ACRES REU student Susannah Yeh
      * Script querys the NCBI non-redundant database through MSU's HPC 
      * Requirements: 
        * BLAST
      * Inputs:
         * `final_prot.fasta`: Protein sequence outpus from xander `search`
         * non-redundant database
      * Outputs:
         * `ncbi.input.txt`: List of accession numbers of top blast hits; used for iTOL tree building
         * `gene.descriptor.final.txt`: List and counts of unique blast descriptor hits; used to see that all gene of interest & minimal hypothetical proteins/ nonspecific targets
         * `accno.final.txt`: Lists accession numbers and their occurrence (ie how many contigs hit to each)
         
## 4. Phylogenetic assessment
  * This step accomplishes the following goals
    * Obtain maximum likelihood trees for each gene of interest
      * One of low quality with no reference sequences (for optional ecological downstream analyses in R)
      * One of higher quality for phylogenetic assessments
    * Construct a gene abundance table for abundance and diversity analyses
  * Calls on [phylo.hmmalign.gene.sh](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/blob/master/assembly_assessments/bin/phylo.hmmalign.gene.sh) unless you want to compile multiple genes into one tree, which calls on [phylo.muscle.group.sh](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/blob/master/assembly_assessments/bin/phylo.muscle.group.sh)
    * Requirements
      * `Bioperl/1.6.923`
      * `raxml/8.0.6`
      * `FastTree/2.1.7`
      * `HMMER/3.1b2`
    * Inputs: 
      * `final_prot_aligned.fasta`: Aligned protein sequence outputs from xander `search`
      * `coverage.txt`: Coverage information for each assembled protein ouptus from xander `search`
      * `get_OTUabundance.sh`: Script to dereplicate and cluster xander output sequences; written by Qiong Wang from the RDP
      * `.seeds`: Seed sequences downloaded from FunGene and used to make hmms
      * `root`: Fasta file containing sequences to root specific tree
      * `matchlist.fasta`: output from BLAST against nr
     * Outputs: 
      * `RAxML_bipartitionsBranchLabels.${GENE}_${CLUST}_RAxML_PROTGAMMAWAG.nwk`: RAxML newick file (100 boostraps)
      * `${GENE}_${CLUST}_FastTree.nwk`: FastTree newick file 
      * `complete.clust_rep_seqs_${CLUST}_unaligned.fasta`: list of representative sequences (at particular cluster)

      

