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

# Troubleshooting
* Cen13 came back with the following error message: 
```
[    844.6]  WARNING: The curve reached near-saturation, hence coverage estimations could be unreliable
 [    844.6]  To avoid saturation increase the -L parameter, currently set at 50.000000
 [    844.6]  WARNING: The curve reached near-saturation in 6 or less points, hence diversity estimations could be unreliable
 [    844.6]  To increase the resolution of the curve increase the -i parameter, currently set at 0.010000
 ```
    * I believe this is because I thought reads were already paired and split them? I will try on the "unsplit" sample
    * If this does not adjust the error message, I will change the suggested parameters
 
* Cen06
  * Sample did not give initial error messages, but threw a warning when plotting saying that the slope is too high; error message suggested decreasing i parameter
  * Will change i from `0.01` to `0.001`
  * Sample is being re-run `nonpareil-mpi -s Cen01.1.fasta -f fasta -t 8 -i 0.001 -b Cen01.txt -R 4194303`
