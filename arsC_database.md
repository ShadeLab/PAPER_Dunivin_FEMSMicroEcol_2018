
#arsC database construction
### Taylor Dunivin
## January 20, 2016
### Goals: 
* Determine distinguishing features between types of *arsC*
* Determine quality of FunGene database
* Distinguish *arsC* from *spxA*

### *arsC* distinguishing features
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

### Current arsC FunGene database
| Accession | Name | Prot Family | Lenth (aa) |
| --------- | ----- | ------ | :-----: |
| Q9CGN5_LACLA (Q9CGN5) | ykhE hypothetical protein | arsC_spx | 115 |
| YQGZ_BACSU (P54503) | |
| Q986E7_RHILO (Q986E7) | |
| Q9A7S8_CAUCR (Q9A7S8) | |
| Q9HXX5_PSEAE (Q9HXX5) | |
| Y103_HAEIN (P44515) | |
| Q9CM21_PASMU (Q9CM21) | |
| YFFB_ECOLI (P24178) | |
| Q9KQ54_VIBCH (Q9KQ54) | |
| Q9JZ90_NEIMB (Q9JZ90) | |
| Q9PH31_XYLFA (Q9PH31) | |
| YUSI_BACSU (O32175) | |
| Q9K785_BACHD (Q9K785) | |
| Q9RY16_DEIRA (Q9RY16) | |
| Q9A861_CAUCR (Q9A861) | |
| Q9CJQ0_PASMU (Q9CJQ0) | |
| Q9A083_STRP1 (Q9A083) | |
| Q98J03_RHILO (Q98J03) | |
| ARSC_SHIFL (P0AB97) | |
| Q9CLS9_PASMU (Q9CLS9) | |
| Q9KQ39_VIBCH (Q9KQ39) | |
| YFGD_ECOLI (P76569) | |
| Q9I508_PSEAE (Q9I508) | |
| ARSC_NEIMB (P63622) | |
| SPX_BACSU (O31602) | |
| SPX_BACHD (Q9K8Z1) | |
| SPX_STRP1 (P60381) | |
| SPX1_LACLA (Q9CI20) | |
| SPX_LACLM (P60376) | |
| Q99XP2_STRP1 (Q99XP2) | |
| Q9CFY7_LACLA (Q9CFY7) | |
| Q9CH93_LACLA (Q9CH93) | |
| Q9CFZ8_LACLA (Q9CFZ8) | |
| Y176_UREPA (Q9PQW7) | |
| Y127_MYCGE (P47373) | |
| Y266_MYCPN (P75509) | |
