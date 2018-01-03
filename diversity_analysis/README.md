# Analyzing diversity of genes *across* chronosequence
### Taylor Dunivin

## Data analysis using R 
   * [post_assembly.R](https://github.com/ShadeLab/PAPER_Dunivin_Antibiotics_2017/blob/master/diversity_analysis/post_assembly.R) contains all R workflows  
      * Calls on several outputs
         * xander `search` outputs
         * BLAST results
         * gene abundance tables (gene-based OTU tables)
         * taxanomic abundance (for rplB)
      * Also includes some manually currated files
         * Site metadata
         * Spatial distance matrix (spacial distance between sites from previous publication [Lee and Sorensen 2016](https://www.ncbi.nlm.nih.gov/pubmed/28282042))
         * Gene catergories (with colors for figure consistency)
         * Phylum colors 
   * All inputs are included in /data
   * Non-figure outputs are included in /outputs
   * Figures are included in /figures
   * Goals of the analysis
      * What is the abundance of assembled ARGs along the Centralia chronosequence?
      * How does diversity change?
  
