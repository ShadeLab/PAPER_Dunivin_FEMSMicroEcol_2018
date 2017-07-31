# Notes on blasting clustered OTUs against JGI assembled contigs

```
module load BLAST+

makeblastdb -in /mnt/scratch/dunivint/Cen06.scaffolds.fasta  -dbtype nucl -out /mnt/scratch/dunivint/databases/Cen06.db

blastn -db /mnt/scratch/dunivint/databases/Cen06.db -query clustered_AsRG_ARG_nucl.fa -out Cen06_blast.txt -outfmt "6 qseqid salltitles pident evalue" -max_target_seqs 10
```

# Notes on blasting clustered OTUs against MAGs

```
cat /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning/VerySpecific_Bins/* > METABAT.fa

module load BLAST+

makeblastdb -in /mnt/scratch/dunivint/METABAT.fa  -dbtype nucl -out /mnt/scratch/dunivint/databases/METABAT.db

blastn -db /mnt/scratch/dunivint/databases/METABAT.db -query /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/networks/clustered_AsRG_ARG_nucl.fsa -out metabat_blast.txt -outfmt "6 qseqid salltitles pident evalue" -max_target_seqs 10
```

Once top hits are generated, you can search for the specific contig in the original files
```
cd /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning/VerySpecific_Bins
find . -name "CONTIG_NAME" -type f -delete
```

This list will let you compare top hits to see if any occur in the same bin. 

| GENE | OTU | Contig | % identity | e-value | Bin | 
| ---- | --- | ------ | ---------- | ------- | --- |
|acr3 | OTU_0053 | k107_5410911 |	98.73	| 0.0 | 608 |
|acr3 | OTU_0053	| k107_374618	| 88.89	| 0.0 | 608, 795|
|acr3 | OTU_0053	| k107_30885960	| 86.16	| 0.0 | 795 |
|acr3 | OTU_0053	| k107_3117242	| 84.88	| 0.0 | 559 |
|acr3 | OTU_0053	| k107_24027065	| 81.54	| 0.0 | 736 |
|acr3 | OTU_0053	| k107_2881296	| 77.52	| 1e-164 | 251|
|acr3 | OTU_0053	| k107_20985308	| 76.86	| 6e-153 | 795 |
|acr3 | OTU_0053	| k107_6209637	| 73.29	| 2e-92 | 666 |
|acr3 | OTU_0053	| k107_17335913	| 71.85	| 7e-58 | 231 |
|acr3 | OTU_0053	| k107_22134734	| 73.06	| 9e-57 | 897 |
|acr3 | OTU_0002	| k107_14595991	| 89.19	| 0.0 | 419 |
|acr3 | OTU_0002	| k107_14536628	| 88.13	| 1e-171 | 827 |
|acr3 | OTU_0002	| k107_17642053	| 73.58	| 3e-14 | 654 |
|arsM | OTU_0109	| k107_1005213	| 84.23	| 7e-135 | 632 |
|arsM |OTU_0109	| k107_35201236	| 77.38	| 3e-78 | 848 |
|arsM |OTU_0109	| k107_3530183	| 76.36	| 4e-38 | 742 |
|arsM |OTU_0109	| k107_33012804	| 76.64	| 1e-27 | 783 |
|arsM |OTU_0296	| k107_19011909	| 77.34	| 4e-119 | 1040 |
|arsM |OTU_0296	| k107_10364269	| 78.57	| 7e-112 | 446 |
|arsM |OTU_0296	| k107_1728533	| 78.48	| 2e-107 | 215 |
|arsM |OTU_0296	| k107_9939490	| 78.40	| 7e-107 | 319 |
|arsM |OTU_0296	| k107_9261249	| 78.93	| 4e-99 | 45 |
|arsM |OTU_0296	| k107_33012804	| 79.43	| 3e-90 | 783 |
|arsM |OTU_0296	| k107_1593047	| 76.92	| 7e-87 | 848 |
|arsM |OTU_0296	| k107_1871211	| 76.21	| 7e-87 | 27 |
|arsM |OTU_0296	| k107_17033332	| 75.55	| 2e-78 | 790 |
|arsM |OTU_0296	| k107_13293973	| 75.59	| 5e-78 | 559|
|tolC |OTU_0005	| k107_11533305	| 76.50	| 3e-156 | 232 |
|tolC |OTU_0005	| k107_15504113	| 79.00	| 4e-155 | 1040 |
|tolC |OTU_0005	| k107_32829128	| 76.02	| 2e-117 | 754 |
|tolC | OTU_0005	| k107_10071582	| 74.83	| 7e-98 | 479 |
|ClassB | OTU_0003	| k107_2170412	| 91.45	| 0.0 | 795 |
|ClassB | OTU_0003	| k107_15390298	| 82.17	| 0.0 | 527 |
|drfa12 | OTU_0076	| k107_16275153	| 98.55	| 1e-62 | 966 |
|drfa12 |OTU_0076	| k107_2833959	| 86.13	| 8e-34 | 591 |
|drfa12 |OTU_0076	| k107_15496897	| 84.06	| 2e-29 | 591 |
|drfa12 |OTU_0076	| k107_11473637	| 84.44	| 2e-29 | 591 |
|drfa12 |OTU_0076	| k107_24297536	| 85.71	| 2e-29 | 361 |
|drfa12 |OTU_0076	| k107_2230677	| 85.04	| 3e-28 | 591|
|drfa12 |OTU_0076	| k107_35661732	| 78.83	| 5e-16 | 298|
