# arsB database construction
### Taylor Dunivin
## February 8, 2017
### Goals: 
* Determine quality of _arsB_ FunGene database
* Determine distinguishing features between *arsB* and *acr3*
* Begin listing seed sequences for *arsB*

### FunGene database quality
* My first thoughts on the db are that there are very few sequences. I think there is more quality diversity to include in the HMM for *arsB* today
* The table below shows each existing seed sequence, its NCBI hit name, its NCBI groups (multidomain, protein superfamily, specific group), length of sequence, and the UniProt annotation score (1-5). 

| Accession | Name | Multidomain | Super Family | Specific group | Lenth (aa) | UniProt Score |
| --------- | ----- | ---------- | --------- | -------- | ------- | :-----: |
| ARSB_STAXY (Q01255) | Arsenical pump membrane protein | arsB | ArsB_NhaD permease | PRK15445 | 429 | 2 |
| YDFA_BACSU (P96678) | Arsenical pump membrane protein | arsB | ArsB_NhaD permease | | 435 | 2 |
| O68021_PSEAE (O68021) | Arsenical pump membrane protein | arsB | ArsB_NhaD permease |  | 425 | 2 |
| ARSB2_ECOLX (P52146) | Arsenical pump membrane protein | arsB | ArsB_NhaD permease | | 429 | 2 |

### *arsB* Notes from FunGene searching
* all existing sequences have low scores and are inferred only by protein homology
* rather than catalytic residues, distinguishing features for *arsB* will likely be transmembrane domains (number and location)
* will need to distinguish between NhaD and ArsB

