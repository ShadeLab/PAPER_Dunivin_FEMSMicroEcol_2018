## Github Repository for
# Community structure explains antibiotic resistance gene dynamics over a temperature gradient in soil
## by Taylor K Dunivin and Ashley Shade


<i>This work is published.</i>


### Data
The raw data for this study are available from [IMG/GOLD study ID: Gs0114513](https://gold.jgi.doe.gov/biosamples?Study.GOLD%20Study%20ID=Gs0114513)


### To cite this work or code
Dunivin TK and A Shade. 2018. Community structure explains antibiotic resistance gene dynamics over a temperature gradient in soil. FEMS Microbiology Ecology (Environmental Dimensions of Antibiotic Resistance special issue). 94: https://doi.org/10.1093/femsec/fiy016



### Abstract
Soils are reservoirs of antibiotic resistance genes (ARGs), but environmental dynamics of ARGs are largely unknown. Long-term disturbances offer opportunities to examine microbiome responses at scales relevant for both ecological and evolutionary processes and can be insightful for studying ARGs. We examined ARGs in soils overlying the underground coal seam fire in Centralia, PA, which has been burning since 1962. As the fire progresses, previously hot soils can recover to ambient temperatures, which creates a gradient of fire impact. We examined metagenomes from surface soils along this gradient to examine ARGs using a gene-targeted assembler. We targeted 35 clinically relevant ARGs and two horizontal gene transfer-related genes (intI and repA). We detected 17 ARGs in Centralia: AAC6-Ia, adeB, bla_A, bla_B, bla_C, cmlA, dfra12, intI, sul2, tetA, tetW, tetX, tolC, vanA, vanH, vanX and vanZ. The diversity and abundance of bla_A, bla_B, dfra12 and tolC decreased with soil temperature, and changes in ARGs were largely explained by changes in community structure. We observed sequence-specific biogeography along the temperature gradient and observed compositional shifts in bla_A, dfra12 and intI. These results suggest that increased temperatures can reduce soil ARGs but that this is largely due to a concomitant reduction in community-level diversity.

### Funding
Metagenome sequencing was supported by the Joint Genome Institute Community Science Project #1834. The work conducted by the U.S. Department of Energy (DOE) Joint Genome Institute, a DOE Office of Science User Facility, is supported under Contract No. DE-AC02-05CH11231. TKD was supported by the Department of Microbiology and Molecular Genetics Ronald and Sharon Rogowski Fellowship for Food Safety and Toxicology Graduate Fellowship.This work was supported in part by Michigan State University through computational resources provided by the [Institute for Cyber-Enabled Research](https://icer.msu.edu/).


### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)


# Antibiotics_Centralia
Repository contains data and workflows for gene-targeted assembly of antibiotic ressitance genes. 

## 1. Reference gene database construction ([reference_sequences](https://github.com/ShadeLab/PAPER_Dunivin_Antibiotics_2017/tree/master/reference_sequences))
- FunGene databases for ARGs were used to compile diverse, near-full length sequences 
- Search parameters are noted in the [reference_sequences](https://github.com/ShadeLab/PAPER_Dunivin_Antibiotics_2017/tree/master/reference_sequences) repository
- Both protein and nucleotide sequences were downloaded

## 2. Gene targeted assembly
- Assembly was performed using [Xander](https://github.com/rdpstaff/Xander_assembler) from the RDP's 2015 [publication](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6)

## 3. [Assembly assessments](https://github.com/ShadeLab/PAPER_Dunivin_Antibiotics_2017/tree/master/assembly_assessments)
- Examine assembly statistics 
- Examine top hits (to reference databases from step 2 as well as NCBI non-redundant database)
- Make maximum likelihood trees (and corresponding labels) for manual assessment
- Output coverage-adjusted gene abundance tables at relevant clustering identities (90, 97, 99%)

## 4. [Diversity analysis](https://github.com/ShadeLab/PAPER_Dunivin_Antibiotics_2017/tree/master/diversity_analysis)
- Read in coverage-adjusted gene abundance tables to R
- Analyze abundance of different gene groups
- Analyze gene diversity (richness, evenness, beta-diversity)
- Examine correlations between genes and environmetnal variables
