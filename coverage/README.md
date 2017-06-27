# Estimating coverage using nonpareil
### Taylor Dunivin

## 1. Prepare files for nonpareil analysis
   * Calls on script [prep.nonpareil.pl](https://github.com/ShadeLab/Xander_arsenic/blob/master/coverage/bin/prep.nonpareil.pl)
   * To execute, `./prep.nonpareil.pl SITE`
   * Script also uses `FastA.split.pl` from [Luis M. Rodriguez](https://github.com/lmrodriguezr/enveomics/blob/master/Scripts/FastA.split.pl)
   * Inputs: `fastq.gz` sequence file (in my case it is already quality filtered)
   * Outputs: `SITE.1.fasta` and `SITE.2.fasta` which contain unziped, split fasta files for nonpareil analysis
   
## 2. Run nonpareil
   * Example `.qsub` file
   ```
#!/bin/bash -login
 
### define resources needed:
### walltime - how long you expect the job to run
#PBS -l walltime=7:00:00:00
 
### nodes:ppn - how many nodes & cores per node (ppn) that you require
#PBS -l nodes=02:ppn=8
 
### mem: amount of memory that the job will need
#PBS -l mem=300gb
 
### you can give your job a name for easier identification
#PBS -N nonpareil.cen06.i
 
### load necessary modules, e.g.
module load nonpareil
 
### change to the working directory where your code is located
cd /mnt/scratch/dunivint/coverage
 
### call your executable
nonpareil-mpi -s Cen06.1.fasta -f fasta -t 8 -b Cen06.i.txt -R 4194303
```
  * Input: `SITE.1.fasta`
  * Outputs: 
    * `SITE.txt.npa`
    * `SITE.txt.npc`
    * `SITE.txt.npl`
    * `SITE.txt.npo`
    
## 3. View nonpareil curves in R
  * Calls on [nonpareil.R](https://github.com/ShadeLab/Xander_arsenic/blob/master/coverage/bin/nonpareil.R)
  * Inputs: 
      * All `SITE.txt.npo` of interest
      * `descriptors_full.txt`: Tab delimited file containing file names, desired names, and line colors (with RBG scale)
  * Output: nonpareil curve!
  
