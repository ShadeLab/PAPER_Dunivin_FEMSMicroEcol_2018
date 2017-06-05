# Assessment by gene
### Taylor Dunivin

__Goals:__
* Be certain that when considering nr database, the gene of interest is still the top hit
* Add to confidence in assembly assessments
* Determine whether/ how to troubleshoot low quality results from specific genes

### Table of contents
* [Setting up workflow](https://github.com/ShadeLab/Xander_arsenic/blob/master/assessment_by_gene.md#setting-up-workflow)
* [blast.qsub](https://github.com/ShadeLab/Xander_arsenic/blob/master/assessment_by_gene.md#blast.qsub)
* [Gene analysis notes](https://github.com/ShadeLab/Xander_arsenic/blob/master/assessment_by_gene.md#gene-analysis-notes)

### __Setting up workflow:__
1. Access HPCC's nr database
```
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB
```

2. Merge all protein fasta files to one `final_prot.fasta`
```
cat *_final_prot.fasta >final_prot.fasta
```

3. Blast nr database: use protein blast `-p blastp`, limit outputs to the two top hits `-b 2 -v 2`, e-values greater than E-6 `e 1e-6`, run with 8 threads `-a 8`
```
blastall -d nr -i final_prot.fasta -p blastp -o blast.txt -b 2 -v 2 -e 1e-6 -a 8
```

4. Make summary file that only shows sequence descriptors (typically gene name)
```
grep '^>' blast.txt > descriptor.blast.txt
```

### __blast.qsub file__
```
#!/bin/bash -login
 
### define resources needed:
### walltime - how long you expect the job to run
#PBS -l walltime=04:00:00
 
### nodes:ppn - how many nodes & cores per node (ppn) that you require
#PBS -l nodes=1:ppn=8
 
### mem: amount of memory that the job will need
#PBS -l mem=50gb
 
### you can give your job a name for easier identification
#PBS -N blast.GENE
 
### load necessary modules, e.g.
module load BLAST/2.2.26
 
### change to the working directory where your code is located
cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/blast/GENE
 
### call your executable
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB
cat *_final_prot.fasta >final_prot.fasta
blastall -d nr -i final_prot.fasta -p blastp -o blast.txt -b 2 -v 2 -e 1e-6 -a 8
grep '^>' blast.txt > descriptor.blast.txt
```

### Gene analysis notes
