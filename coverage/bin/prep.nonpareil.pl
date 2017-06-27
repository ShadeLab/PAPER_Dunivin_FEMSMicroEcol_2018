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
