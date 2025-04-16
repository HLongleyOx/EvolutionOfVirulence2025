#This script solves the quasispecies equation numerically for each moment in time to determine the frequency of each mutation type

#Install libraries
library(tidyverse)
library(deSolve)
library(abind)
library(pracma)

#Set seed
set.seed(1004)

#Load functions
source("tools/many_mutation_model_calc.R")
source("tools/quasispecies_model_full.R")
source("tools/hill_functions.R")

args <- commandArgs(trailingOnly=TRUE) #batch processing should be used if ppossible
# commands are the number 1) number of virus types (m), 2) selection cost defined in within_host_equilibrium.R 3) The number of mutations carried by initial virus type

#set fixed parameter values
max_dur_years <- (0.25 + hill_down_function(10^2) + 0.75) #max duration in yrs 
max_dur <- max_dur_years*365

num_strains <- as.numeric(args[1])
mu <- 3*10^-5 #mutation rate
s <- as.numeric(args[2])
start_i <- as.numeric(args[3])

fullmodel <- quasispecies_model_full(start_i =start_i, m=num_strains,M=num_strains+1,s=s, max_dur=max_dur, dt=1, mu=mu)

file = paste0("withinhostmodel_",num_strains,"_",start_i,".csv") #save to desired folder 

write_csv(fullmodel, file)

