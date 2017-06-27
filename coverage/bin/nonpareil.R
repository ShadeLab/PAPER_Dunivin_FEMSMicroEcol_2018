source('/Users/dunivint/Documents/GitHubRepos/nonpareil/utils/Nonpareil.R')

Nonpareil.curve('/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/coverage/Cen01.txt.npo')
Nonpareil.legend("bottomright")

tools::Rd2txt(tools::parse_Rd('/Users/dunivint/Documents/GitHubRepos/nonpareil/utils/nonpareil/man/Nonpareil.curve.Rd'))
getwd()
setwd("/Users/dunivint/Documents/GitHubRepos/Xander_arsenic/coverage/")
samples <- read.table('descriptors_mini.txt', sep="\t", h=T);
attach(samples);
np <- Nonpareil.curve.batch(File, r=R, g=G, b=B, libnames=Name, modelOnly=TRUE);
Nonpareil.legend('bottomright');
detach(samples)
