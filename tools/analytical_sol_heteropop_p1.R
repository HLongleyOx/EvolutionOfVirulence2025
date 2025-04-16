
#analytical solution to the between host model when we assume host hetrogeneity 
#host_type = number of host types 
#num strains = number of infection types
#model = within host model output
#spvl = spvl of virus type
#E = host effect size
#returns dataframe of of total prevalence, prevalence by mutation type, spvl distribution 


#host_effect = seq(-E, E, length.out = host_types)

analytical_sol_heteropop_p1 <- function(host_types, num_strains, model,spvl,E) {
    FOI_HH = force_of_infection_calculation(
    time = 7805,
    spvl = spvl,
    muts_freq = model,
    num_strains = num_strains + 1,
    E=E)
    return(list(FOI_HH[[2]],FOI_HH[[3]],FOI_HH[[4]], FOI_HH[[6]]))
}
