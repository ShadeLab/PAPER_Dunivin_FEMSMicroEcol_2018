# Notes on running nonpareil on Centralia samples
### Taylor Dunivin

# Script to prepare sequences for nonpareil
* Located in `/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/coverage/bin/prep.nonpareil.pl
* Script completes the following:
  1. Copies JGI metagenome to my scratch space
  2. Unzips metagenome file
  3. Converts from .fastq to .fasta
  4. Splits paired end reads into two files `${SITE}.1.fa` and `${SITE}.2.fa`
* To excecute, run `./prep.nonpareil.pl <SITE>
* Script also uses `FastA.split.pl` from [Luis M. Rodriguez](https://github.com/lmrodriguezr/enveomics/blob/master/Scripts/FastA.split.pl)

```
#!/bin/bash

#you must specify the site name
SITE=$1

#change to scratch directory
cd /mnt/scratch/dunivint/coverage

#Move files to scratch
scp /mnt/research/ShadeLab/Sorensen/JGI_Metagenomes/${SITE}.anqdp.fastq.gz ${SITE}.anqdp.fastq.gz

#unzip file
gunzip ${SITE}.anqdp.fastq.gz 

#Convert to fasta from fastq
cat ${SITE}.anqdp.fastq | paste - - - - | awk 'BEGIN{FS="\t"}{print ">"substr($1,2)"\n"$2}' > ${SITE}.anqdp.fasta

#Split paired end reads
module load perl
/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/coverage/bin/./FastA.split.pl ${SITE}.anqdp.fasta ${SITE} 2
```

# Running nonpareil
* Run the following command on the `<Site>.1.fa` output file from `prep.nonpareil.pl`
* `nonpareil-mpi -s Cen01.1.fasta -f fasta -t 8 -b Cen01.txt -R 4194303`
