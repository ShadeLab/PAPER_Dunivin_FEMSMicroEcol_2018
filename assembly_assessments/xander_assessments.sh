#!/bin/bash

#you have to type the gene and sample names after calling this script, for example "$./testAutomation.sh arsB cen10"
#$1 means the first argument (arsB in the example, so the variable GENE would be arsB), $2 is site name
GENE=$1
SAMPLE=$2

#load blast
module load GNU/4.4.5
module load Gblastn/2.28

#switch into correct directory. (Not sure which directory you want)
cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/assessments
mkdir ${GENE}
cd ${GENE}

#blast xander results against db
#tabular format, show seq id, query id (and description), e-value, only show 1 match
blastn -db /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/assessments/blast_databases/${GENE}/${GENE}_database -query /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/${SAMPLE}/k45/${GENE}/cluster/*_final_nucl.fasta -out /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/${SAMPLE}/k45/${GENE}/cluster/blast.txt -outfmt "6 qseqid salltitles evalue" -max_target_seqs 1

#make a list of reads from *match_reads.fa
grep "^>" /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/${SAMPLE}/k45/${GENE}/cluster/*match_reads.fa | sed '0~1s/^.\{1\}//g' > /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/${SAMPLE}/k45/${GENE}/cluster/matchreadlist.txt

##Extract info for stats calculations in R
#take the framebot stats
grep "STATS" /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/${SAMPLE}/k45/${GENE}/cluster/*_framebot.txt > /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/${SAMPLE}/k45/${GENE}/cluster/framebotstats.txt

# USE R TO CREATE STATISTIC FILES

#start in cluster directory from xander output!
cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/${SAMPLE}/k45/${GENE}/cluster/

# swap and load the version of R you want to use here with module commands
module load GNU/4.9
module load OpenMPI/1.10.0
module load R/3.3.0
 
# Run R Command with input script myRprogram.R
#R script should read in files from the cluster and create 4 new files
R < /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/assessments/bin/assembly_assessmentR.R --no-save

#move R files to databases_${SAMPLE}

mv readssummary2.txt ${SAMPLE}_${GENE}_readssummary.txt
mv ${SAMPLE}_${GENE}_readssummary.txt /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/assessments/${GENE}

mv kmerabundancedist2.png ${SAMPLE}_${GENE}_kmerabundancedist.png
mv ${SAMPLE}_${GENE}_kmerabundancedist.png /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/assessments/${GENE}

mv stats2.txt ${SAMPLE}_${GENE}_stats.txt
mv ${SAMPLE}_${GENE}_stats.txt /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/assessments/${GENE}

mv e.values2.txt ${SAMPLE}_${GENE}_e.values.txt
mv ${SAMPLE}_${GENE}_e.values.txt /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/assessments/${GENE}


#GET GC COUNT OF THIS SAMPLE!

cd /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/assessments

#load perl
module load GNU/4.4.5
module load perl/5.24.1
#get GC count and put output into ${GENE}_gc folder
bin/./get_gc_counts.pl /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/analysis/${SAMPLE}/k45/${GENE}/cluster/${SAMPLE}_${GENE}_45_final_nucl.fasta

mv gc_out.txt ${SAMPLE}_${GENE}_gc_out.txt
mv ${SAMPLE}_${GENE}_gc_out.txt /mnt/research/ShadeLab/WorkingSpace/Dunivin/xander/assessments/gc
