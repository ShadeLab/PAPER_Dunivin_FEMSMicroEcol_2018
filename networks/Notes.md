# Notes on network analysis
### Taylor Dunivin

### Table of contents:
* [June 26, 2017]()



### June 26, 2017
  * Deciding which network method to use
    * Can't use MIC (under 25 samples)
    * Frontiers article recomends and Tiedje ISME uses simple pearson's R for calculation
        * Performed this analysis
        * Caveat: would need >7 samples to construct network that does not require approximation, but this eliminates most OTU's at 90% identity
    * MENA is a good option
        * Made for microarray data (ie can accomodate questions of per genome interest)
        * non-arbitrary threshold
        
  * MENA workflow
    * Used .R script to export dataset
      1) counts normalized to total xander-assembled rplB count
      2) zeros are left blank
      3) contains all tested AsRG and ARG
      4) contains all 13 sample sites
      5) OTUs are clustered at 90% identity 
    * Constricuted network using OTUs present in at least 4 samples
      * 57% sparse
      * 177/416 total (calculated in R)
    * Selected RMT cutoff of 0.310 as it was marked with p value > 0.001
    ```
The number of nodes: 215
The number of links: 5014
The average path: 1.782
R square of power-law: 0.536
```
    * Examined global network properties
```
Network Indexes	blankz(0.310)
Total nodes	215
Total links	5014
R square of power-law	0.536
Average degree (avgK)	46.642
Average clustering coefficient (avgCC)	0.701
Average path distance (GD)	1.782
Geodesic efficiency (E)	0.609
Harmonic geodesic distance (HD)	1.642
Maximal degree	213
Nodes with max degree	rplB_0186
Centralization of degree (CD)	0.785
Maximal betweenness	1022.554
Nodes with max betweenness	rplB_0186
Centralization of betweenness (CB)	0.041
Maximal stress centrality	17844
Nodes with max stress centrality	rplB_0186
Centralization of stress centrality (CS)	0.716
Maximal eigenvector centrality	0.173
Nodes with max eigenvector centrality	rplB_0186
Centralization of eigenvector centrality (CE)	0.114
Density (D)	0.218
Reciprocity	1
Transitivity (Trans)	0.306
Connectedness (Con)	1
Efficiency	0.786
Hierarchy	0
Lubness	1
```
        
