# Analyzing diversity of genes *across* chronosequence
### Taylor Dunivin

### Table of contents
* [Data preparation](https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis#data-preparation)
* [Data analysis using R](https://github.com/ShadeLab/Xander_arsenic/tree/master/diversity_analysis#data-analysis-using-r)

## Data preparation
### Taxanomic abundance
 * Input: ```${SITE}_${GENE}_${KMER_SIZE}_taxonabund.txt```, `Centralia_full_map.txt`, `gene_classification.txt`
 * Remove `lineage match name` and below since it is essentially a duplicate of above section for non-RDP created gene databases
 * Table includes taxon name, abundance, and relative abundance (of total for gene) for gene of interest
 * Product of ` "search" ` in xander

### Prepare OTU table for each gene of interest
* This is completed in the [phylogenetic analysis](https://github.com/ShadeLab/Xander_arsenic/tree/master/phylogenetic_analysis)
* Within phylogenetic analysis script (`phylo.sh`), there is a nested script created by the RDP called `get_OTUabundances.sh` which outputs an R formatted distance OTU table (`rformat_dist.txt`)
* Step can be completed for any clustering distence (need to specify)

## Data analysis using R 
### 1. Diversity analysis
* Calculate alpha and beta diversity post-rarefication for each individual gene of interest
* Inputs: `rformat_dist.txt`, `Centralia_full_map.txt`
* [diversity_analysis_slim.R](https://github.com/ShadeLab/Xander_arsenic/blob/master/diversity_analysis/diversity_analysis_slim.R)


### 2. Taxanomic abundance analysis (match level; xander default)
* Use one of two scripts that differ based on their normalization
    * [meta_arg_rplB.R](https://github.com/ShadeLab/Xander_arsenic/blob/master/diversity_analysis/meta_arg_rplB.R)
    * [meta_arg_census.R](https://github.com/ShadeLab/Xander_arsenic/blob/master/diversity_analysis/meta_arg_census.R)

### 3. Network analysis & OTU level taxanomic abundance analysis
* Inputs: `rformat_dist.txt`, `Centralia_full_map.txt`, `gene_classification.txt`
* [network_analysis.R](https://github.com/ShadeLab/Xander_arsenic/blob/master/diversity_analysis/network_analysis.R)
* Outputs: `otu_table.txt`, `corr_table.rplB.txt`

