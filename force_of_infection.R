
## This script takes the frequencies of every virus type over time and determines the force of infection at any given moment of the infections, which describes the probabilities of which virus type will be transmitted. 
## Also output the between host (bh) equilibrium distribution at the endemmic equilibria 

#Install libraries
library(tidyverse)
library(deSolve)
library(abind)
library(pracma)

#Set seed
set.seed(1004)

#Load functions
source("tools/force_of_infection_calculation.R")
source("tools/viralload_over_time.R")
source("tools/hill_functions.R")
source("tools/bh_equilbria.R")


args <- commandArgs(trailingOnly=TRUE) #preferably batch processing should be used due to size of data 

#set fixed parameter values
max_dur_years <- (0.25 + hill_down_function(10^2) + 0.75) #max duration in yrs 
max_dur <- max_dur_years*365

num_strains <- as.numeric(args[1])
n=num_strains
muts_freq = array(dim=c(n+1,n+1,max_dur))

for (i in 0:n){
  df = read_csv(paste0("withinhostmodel_",n,"_",i,".csv")) #read in frequencies 
  df = df[,c(2:(n+2))]
  for (t in 1:max_dur){
    muts_freq[,(i+1),t] = as.numeric(df[t,]) #reformats matrix
  }
}

save(muts_freq, file=paste0("muts_freq_",n,".RData")) #Write R object

vl_cost <- 5/n
vl <-  7-(0:(n))*vl_cost

vl_time <- viralload_over_time(n, muts_freq, vl)
FOI <- force_of_infection_calculation(time=max_dur, spvl=vl, muts_freq=muts_freq, E=0, n+1 )
durations <- FOI[[2]]
spvl_overTime <- FOI[[3]]
spvl = c()

for (i in 1:(num_strains+1)){
  d = (durations[i]*365)-0.75*365
  spvl[i] = mean(spvl_overTime[(0.25*365):d,i])
}

equil_prev <- bh_equilbria(FOI, durations,  spvl,num_strains+1, n) #bh equilibrium distribution 
filepath=paste0("bh_prev_",n,".RData")
save(equil_prev, file=filepath)


#write FOI to a file 
filepath=paste0("FOI_chronic_",n,".RData")
save(FOI, file=filepath)

#write viral load over time to a file 
write_csv(vl_time, paste0("VL_",n,".csv"))




