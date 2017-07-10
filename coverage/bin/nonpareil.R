#note: script is written for personal computer due to installation specifications

#use nonpareil 
source('/Users/dunivint/Documents/GitHubRepos/nonpareil/utils/Nonpareil.R')

#setwd
setwd("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/coverage/data")

#read in sample information (filename, sample name, colors) and then plot
samples <- read.table('descriptors_full.txt', sep="\t", h=T);
attach(samples);
np <- Nonpareil.curve.batch(File, r=R, g=G, b=B, libnames=Name, modelOnly=TRUE);
detach(samples)
