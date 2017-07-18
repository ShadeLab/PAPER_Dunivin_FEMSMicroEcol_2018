# Analyzing diversity of genes *across* chronosequence
### Taylor Dunivin

## Data preparation
### Taxanomic abundance
 * Use ```${SITE}_${GENE}_${KMER_SIZE}_taxonabund.txt``` 
 * Remove `lineage match name` and below since it is essentially a duplicate of above section for non-RDP created gene databases
 * Table includes taxon name, abundance, and relative abundance (of total for gene) for gene of interest
 * Product of ` "search" ` in xander

### Prepare OTU table for each gene of interest
* This is completed in the [phylogenetic analysis](https://github.com/ShadeLab/Xander_arsenic/tree/master/phylogenetic_analysis)
* Within phylogenetic analysis script (`phylo.sh`), there is a nested script created by the RDP called `get_OTUabundances.sh` which outputs an R formatted distance OTU table (`rformat_dist.txt`)
* Step can be completed for any clustering distence (need to specify)

## Data analysis using R 
### 1. Diversity analysis (alpha and beta)


### 2. Taxanomic abundance analysis (match level; xander default)


### 3. Network analysis


### 4. OTU level taxanomic abundance analysis
