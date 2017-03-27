# Compiling gene references working protocol
### Taylor Dunivin
## March 27, 2017
---
## Linking protein and nucleotide accession/GI numbers
### Format sequence information (protein)
```
# make new file of accession number only
grep "^>" input.fa | sed '0~1s/^.\{1\}//g'| cut -f1 -d " "  >prot.id.final.txt
```

### Submit job to assign nucleotide information to protein information
Below is the information in the file ```job.names.qsub```

```
#!/bin/bash -login
 
### define resources needed:
### walltime - how long you expect the job to run
#PBS -l walltime=03:00:00
 
### nodes:ppn - how many nodes & cores per node (ppn) that you require
#PBS -l nodes=10:ppn=1
 
### mem: amount of memory that the job will need
#PBS -l mem=5gb
 
### you can give your job a name for easier identification
#PBS -N prot2accession

### change to the working directory where your code is located
cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/intI
 
### call your executable
./fetchCDSbyProteinIDs_ver2.py prot.id.final.txt names.fa names.txt
```

To submit the job, ```qsub names.qsub```

### Dereplicate nucleotide sequences
First need to dereplicate based on accession number (otherwise RDP software will not work)
```
#remove duplicate accession numbers
awk '/^>/{f=!d[$1];d[$1]=1}f' input.fa >derepaccno.input.fa
```

Next dereplicate based on sequence
```
java -Xmx2g -jar /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/RDPTools/Clustering.jar derep -o derep.nucl.fa derep.all_seqs.ids derep.all_seqs.samples derepaccno.input.fa
```

### Obtain accession numbers for nucleotide sequences in ```derepaccno.fa```
```
# make new file of accession number only
grep "^>" derep.nucl.fa | sed '0~1s/^.\{1\}//g'| cut -f1 -d " "  >derep.nucl.id.txt
```

### Remove sequences in protein fasta based on derep nucleotide sequences
Here i will use the filterbyname.sh funciton in bbmap
```
module load BBMap/35.34
filterbyname.sh in=input.fa out=prot.filtered.fa names=derep.nucl.id.txt
```


