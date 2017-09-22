# ARG-AsRG_co-occurrence_Centralia
Repository contains data and workflows for gene-targeted assembly of antibiotic and arsenic resistance genes (ARG and AsRG respectively). 

## 1. AsRG database construction ([seed_sequences](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/tree/master/seed_sequences))
- [seed_sequences](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/tree/master/seed_sequences) were compiled and submitted to the RDP
- Databases are publicly available at [FunGene](http://fungene.cme.msu.edu/)

## 2. Reference gene database construction ([reference_sequences](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/tree/master/reference_sequences))
- FunGene databases for ARG and AsRG were used to compile diverse, near-full length sequences 
- Search parameters are noted in the [reference_sequences](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/tree/master/reference_sequences) repository
- Both protein and nucleotide sequences were downloaded

## 3. Gene targeted assembly
- Assembly was performed using [Xander](https://github.com/rdpstaff/Xander_assembler) from the RDP's 2015 [publication](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6)

## 4. [Assembly assessments](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/tree/master/assembly_assessments)
- Examine assembly statistics 
- Examine top hits (to reference databases from step 2 as well as NCBI non-redundant database)
- Make maximum likelihood trees (and corresponding labels) for manual assessment
- Output coverage-adjusted gene abundance tables at relevant clustering identities (90, 97, 99%)

## 5. [Diversity analysis](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/tree/master/diversity_analysis)
- Read in coverage-adjusted gene abundance tables to R
- Analyze abundance of different gene groups
- Analyze gene diversity (richness, evenness, beta-diversity)
- Examine co-occurrence patterns between genes and environmetnal variables
- Network analysis of ARG and AsRG 

## 6. [Phylogenetic analysis](https://github.com/ShadeLab/ARG-AsRG_co-occurrence_Centralia/tree/master/phylogenetic_analysis)
- Use [iTOL](http://itol.embl.de/) for visualization 
- Add publication-ready labels
- Add abundance (from step 5)

