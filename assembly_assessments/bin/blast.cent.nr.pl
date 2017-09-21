#!/bin/bash

#you have to type the gene and project
GENE=$1
BASEDIR=/mnt/research/ShadeLab/Dunivin/gene_targeted_assembly
RESULTDIR=/mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis

#change to working directory
cd ${BASEDIR}/centralia_arsenic/ncbi_blast
mkdir ${GENE}
cd ${GENE}

#copy all files to database
for i in cen01 cen03 cen04 cen05 cen06 cen07 cen10 cen12 cen13 cen14 cen15 cen16 cen17; do cp ${RESULTDIR}/${i}/k45/${GENE}/cluster/${i}_${GENE}_45_final_prot.fasta ${BASEDIR}/centralia_arsenic/ncbi_blast/${GENE}/${i}_${GENE}_45_final_prot.fasta ; done

#concatenate all fasta files
cd ${BASEDIR}/centralia_arsenic/ncbi_blast/${GENE}
cat *${GENE}_45_final_prot.fasta > final_prot.fasta

#blast nr database!
module load GNU/4.9
module load BLAST+/2.6 
export BLASTDB=/mnt/research/ShadeLab/Dunivin/gene_targeted_assembly/blast_databases/ncbi:$BLASTDB
blastp -db nr -query final_prot.fasta -out blast.txt -outfmt "6 qseqid stitle sacc evalue score length pident" -num_threads 8 -max_target_seqs 1
blastp -db nr -query final_prot.fasta -out gene.descriptor.txt -outfmt "6 stitle" -num_threads 8 -max_target_seqs 1 
blastp -db nr -query final_prot.fasta -out accno.txt -outfmt "6 sacc" -max_target_seqs 1 -num_threads 8

#count occurrences of gene descriptions
sort gene.descriptor.txt | uniq -c > gene.descriptor.final.txt

#get and count occurrences of accession numbers
sort accno.txt | uniq > accno.u.txt

filename='accno.u.txt'
echo Start
while read p; do 
   curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$p&rettype=fasta&retmode=text" >$p.result
done <$filename

#convert results into one matchlist fasta file
cat *.result >matchlist.fasta

#remove unnecessary files
rm *final_prot.fasta
rm *.result
rm accno.txt
rm gene.descriptor.txt
