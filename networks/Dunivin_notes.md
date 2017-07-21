# Notes on blasting clustered OTUs against JGI assembled contigs

```
module load BLAST+

makeblastdb -in /mnt/scratch/dunivint/Cen06.scaffolds.fasta  -dbtype nucl -out /mnt/scratch/dunivint/databases/Cen06.db

blastn -db /mnt/scratch/dunivint/databases/Cen06.db -query clustered_AsRG_ARG_nucl.fa -out Cen06_blast.txt -outfmt "6 qseqid salltitles pident evalue" -max_target_seqs 10
```
