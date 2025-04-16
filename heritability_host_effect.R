
#This script estimates heritability when there are 50 different host types, each of which modifies the viral load to a different extent

library(tidyverse)


heritability_host_effect <- function(num_strains, host_effect, host_types=50){

load(paste0("hetro_bh_", num_strains-1, "_", host_effect, ".RData"))

#create additional dataset that sums over the same infection type but in different hosts 
  prev_spvl <- analypart2[[1]]
  km <- analypart2[[2]]
  sum_inf_mat <- matrix(NA, nrow = num_strains, ncol = num_strains*host_types)

  for (i in 1:num_strains){
    selected_rows <- seq(i, nrow(km), by = num_strains) #identify all type j infections across host types
    sum_inf <- colSums(km[selected_rows,]) 
  
  #add to new matrix
   sum_inf_mat[i,] <- sum_inf
  }


  viral_loads <- data.frame()
  heritability_df <- data.frame()

#loop through 
  for (j in 1:1000){
    for (i in 1:100){
#sample the host type (random number between 1 and 50)
      host_type_source <- sample(1:host_types, size=1)

###sample source infection type based upon the prevalence at equilibrium
#Take only the rows related to the host type 
      selected_rows <- seq((host_type_source-1)*num_strains + 1, (host_type_source-1)*num_strains + num_strains)
      prev_of_host_type <- prev_spvl[selected_rows,]
      source_inf_type <- sample(1:num_strains, size=1, prob=abs(prev_of_host_type$prevalence))
      source_inf_type <- sample(1:num_strains, size=1)
      source_spvl <- prev_of_host_type[source_inf_type,1]
      print(source_spvl)
#based upon the infection of the source, sample the infection of the recipient based upon the transmission potential 
      p <- sum_inf_mat[,(num_strains*(host_type_source-1))+source_inf_type]
      rec_inf_type <- sample(1:num_strains, 1, prob=p)
      rec_host_type <- sample(1:host_types, size=1)
      rec_spvl <- prev_spvl[(num_strains*(rec_host_type-1))+rec_inf_type,1]
      print(c(source_spvl, rec_spvl))
      viral_loads <- rbind(viral_loads,data.frame(rec_spvl=rec_spvl, source_spvl=source_spvl))
  }
    r <- summary(lm(rec_spvl~source_spvl, data=viral_loads))$coeff[[2]]
   # print(summary(lm(rec_spvl~source_spvl, data=viral_loads)))
    heritability_df <- rbind(heritability_df, data.frame(herit_estimate=r, num_strains=num_strains, host_effect=host_effect))
    }

  return(heritability_df)
}


n <- c(11,51,101,151,201,251)
m <- c(0.1,0.25,0.5,0.75,1)
heritability_all <- data.frame()

for (i in n){
  for (j in m){
    
    output <- heritability_host_effect(i, j)
    heritability_all <- rbind(heritability_all, output)
     
  }
}


write_csv(heritability_all, "heritability_df.csv")

