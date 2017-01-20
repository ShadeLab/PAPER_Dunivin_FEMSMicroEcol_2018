
#arsC database construction
### Taylor Dunivin
## January 20, 2016
### Goals: 
* Determine distinguishing features between types of arsC
* Determine quality of FunGene database
* Distinguish arsC from spxA

### arsC distinguishing features
| NCBI protein family       | Type     | Feature 1      |  2 |  3 |  4 |  5 | Length (aa) |  Comments  |
| ------------- | ----- | ----- | ----- | ----- | ----- | ---- | ------ | :---------------------: |
| arsC_15kD    | Glutaredoxin | C (10) | |  N (56) | R (90) | F (103)| 111-112 | COG1393, nitrogenase assoc? |
| arsC_arsC   | Glutaredoxin | C/X (9) | | R (58) | R (92) | R (105) | 111 | crystal incl |
| arsC_YffB   | Glutaredoxin | C (9) | | N (56) | R (92) | F (106) | 104-106 | crystal incl |
| arsC_like  | Glutaredoxin | C (9) | | N (60) | R (96) | F (110) | 109-112 | |
| arsC_spx  | na | C (10) | C/S (13) | | | | 115-120| no confirmed As-relationship|
| arsC_pI258_fam | Thioredoxin | | | | | | 122-128 | crystal incl |

From the above table, I am now skeptical of the arsC_spx family as actually containing arsC since only one catalytic site is shared and there is no experimental evidence from the NCBI sequences presented to date that sequences from this protein family is an arsenate reductase. 

The below tree shows the arsC (glutaredoxin) sequences according to NCBI. Color code is as follows
* Yellow: arsC_arsC
* Grey: arsC_15kD
* Lime: arsC_spx
* Orange: arsC_like
* Blue: arsC_Yffb

![Image of arsC tree](https://github.com/ShadeLab/Xander_arsenic/blob/master/arsC_family_tree.gif)

Even the tree of the arsC (glutaredoxin) suggests that arsC_15kD and arsC_arsC are indeed separate from other protein families. With the above information considered, I would say that only protein families arsC_15kD and arsC_arsC are indeed indicative of glutaredoxin arsenate reductases. 
