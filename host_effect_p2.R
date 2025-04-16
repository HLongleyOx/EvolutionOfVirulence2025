
#Script determines the frequencies of each mutation type when there are different host types which modify the viral load associated with the virus
library(tidyverse)
library(deSolve)
library(abind)
library(pracma)

#FOI

args <- commandArgs(trailingOnly=TRUE) #batch processing if possible. arg 1 is the number of virus types (m) and arg 2 is the maximum size of the host effect on log viral load 

source("tools/analytical_sol_heteropop_p2.R")
source("tools/force_of_infection_calculation.R")
source("tools/hill_functions.R")

#set fixed parameter values
num_strains <- as.numeric(args[1])
host_types <- 50
E <- as.numeric(args[2])
vl_cost <- 5/num_strains
spvl <-  7-seq(0,num_strains)*vl_cost


analypart2  <- analytical_sol_heteropop_p2(host_types, num_strains, spvl,E)
save(analypart2,file=paste0("hetro_bh_", num_strains,"_",E,".RData"))


