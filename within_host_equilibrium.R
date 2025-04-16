#This scripts sets the selection coefficients for our range of interest, where values are chosen so that m=250 has a fitness of 5e-5
#Outputs dataset with equilibrium distribution of virus types 


#Install libraries
library(tidyverse)
library(deSolve)
library(abind)
library(pracma)


#Load functions
source("tools/quasispecies_model_equilibrium_full.R")


## Set selection coefficient 
#fix m~s relationship such that 
mu  = 3 * 10^-5
m <- seq(10, 500, by=10)

get_selection_coefficient <- function(x) {
  # Parameters for line from x=10,y=1e-2 to x=250,y=5e-5
  m = -0.009587625  # slope on log10 scale
  b = -1.904124     # intercept on log10 scale
  
  # Calculate s: log10(s) = mx + b
  s = 10^(m * x + b)
  
  return(s)
}

fullmodel_return=data.frame()    
for (num_strains in m){
s=get_selection_coefficient(num_strains)
fullmodel_equilibrium1 <- quasispecies_model_equilibrium_full(m=num_strains,M=num_strains+1,s=s, mu=mu) %>%
  dplyr::mutate(s=s, m=num_strains)
vl_cost <- 5/num_strains
vl_t_10 <- 7-(0:num_strains)*vl_cost
fullmodel_equilibrium1$vl = sum(vl_t_10*fullmodel_equilibrium1$equilibrium)
fullmodel_return=rbind(fullmodel_return, fullmodel_equilibrium1)
}

write_csv(fullmodel_return,"equilib_vl.csv") #save to desired output folder 
