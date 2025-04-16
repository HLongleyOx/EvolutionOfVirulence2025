
# This script gives Heritability estimate by sampling source infections by using equilibrium prevalence then sampling a corresponding recipient based up Km matrix that describes expected number on onward infections for each virus type

#Sample 1000 infections based upon equilibrium prevalence 

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
num_strains <- as.numeric(args[1])

set.seed(1004)

load(paste0("bh_prev_", num_strains,".RData"))
load(paste0("FOI_chronic_", num_strains,".RData"))

prev = equil_prev[[1]]$prevalence
herit_df <- data.frame()

durations <- FOI[[2]]
spvl_overTime <- FOI[[3]]
spvl = c()

for (i in 1:(num_strains+1)){
  d = durations[i]*365
  spvl[i] = mean(spvl_overTime[1:d,i])
}

herit_coeff <- data.frame()

for (i in 1:100){
  infs <- sample(0:num_strains, 100, replace=T, prob=prev)
  KM =FOI[[4]]

#for each individual, sample a follow on infection based upon the KM 
  for (inf in infs){
    probs <- KM[, (inf+1)] * (1/sum(KM[,(inf+1)]))
    recip <- sample(0:num_strains, 1, prob=probs)
  
    spvl_donor <- spvl[inf+1]
    spvl_recip <- spvl[recip+1]
  
    herit_df <- rbind(herit_df, data.frame(spvl_donor=spvl_donor, spvl_recip=spvl_recip, i=i))
  }
  herit <- summary(lm(spvl_recip~spvl_donor, data=herit_df))$coefficients[[2]]
  herit_coeff <- rbind(herit_coeff, data.frame(herit,i))
}

write_csv(herit_df, file=paste0("heritdf_", num_strains,"_random.csv"))




