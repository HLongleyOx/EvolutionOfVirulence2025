library(deSolve)
library(abind)
library(pracma)

#Set seed
set.seed(1004)

#Load functions
source("tools/between_host_model.R")
source("tools/hill_functions.R")


args <- commandArgs(trailingOnly=TRUE) #should use batch processing due to computational time
#first argument is number of virus types (m) and second is the starting number of viral mutations
#set fixed parameter values
max_dur_years <- (0.25 + hill_down_function(10^2) + 0.75) #max duration in yrs 
max_dur <- max_dur_years*365

num_strains <- as.numeric(args[1])
i <- as.numeric(args[2])
n=num_strains

#Read in force of infection data and virus type frequency 
load(paste0("FOI_",n,".RData"))
load(paste0("muts_freq_",n,".RData"))

durations <- FOI[[2]]
time_taken = 50 #time step
dt=time_taken/365
AICint = max_dur_years  #Duration of initial conditions: the max duration of the infection
iT0 = ceil(AICint/dt) #Time in days
num_host_types = 1
inds <- c(TRUE,rep(FALSE, time_taken-1))
fA = FOI[[1]][,,inds]
mod = list()

#start with high viral load/average viral load/low viral load
initialVL= data.frame(n/5,2*n/5 ,3*n/5, 4*n/5, 5*n/5 )

mod = between_host_model(f=fA,
                                 tmax=200,
                                 tstep=time_taken,
                                 max_dur=max_dur,
                                 initCons=initialVL[,i],
                                 num_strains=n+1,
                                 num_hosts=1,
                                 durations_i=durations,
                                 host_dist=1,
                                 transmissibility=1)


save(mod,file=paste0("mod_",n,"_",i,".RData")) 

