#!/bin/bash

GENE=$1
CLUST=$2
PROJECT=$3
BASEDIR=/mnt/research/ShadeLab/Dunivin/gene_targeted_assembly

cd ${BASEDIR}/${PROJECT}/OTUabundances/${GENE}
grep "^>" complete.clust_full_rep_seqs_${CLUST}_unaligned_seeds.fasta | sed '/^>OTU/ d' | sed 's/[>]//g' | sed -e :1 -e 's/\(\[[^]]*\)[[:space:]]/\1_/g;t1' | awk '{print $1",", $NF,"("$1")"}' | sed 's/\[//g' | sed 's/\]//g' > labels.txt
cat ${BASEDIR}/bin/iTOL.txt labels.txt >labels_iTOL_${CLUST}_${GENE}.txt
