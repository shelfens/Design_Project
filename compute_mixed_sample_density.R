#clear all variables :
  rm(list = ls())

#library for traingular distribution :
  library(triangle)

#EXAMPLE FOR CAMPYLOBACTER ,study of Petterson et al. 2016: 

X <- 500000             # fictive population size 
incidence <- 2.4/100    # [%] (for campylobacter = 2.4%)
rate <- 250             # [g/cap/day] 

# Yearly production of non-infected feces by the non infected population :
prod_no_infected_pop <- 250*365*(1-incidence)*X  # [g/year] 


# Yearly production of infected feces by the infected population :

# Set parameters for Duration [days]
min_duration <- 15   # Minimum value
max_duration <- 42   # Maximum value
mode_duration <- 34  # Mode (most probable value)
n <- X*incidence  # Number of random numbers to generate (= infected population)

# Generate random numbers from a triangular distribution : 
# attribute a duration of excretion [days] to each infected person
random_duration <- rtriangle(n, min_duration, max_duration, mode_duration)
random_duration_rounded = ceiling(random_duration)
  # We need to round to nearest whole number the values of duration
  # (to have entires days)

#Compute the yearly production of infected feces by the infected population : 
  prod_infected_pop_infected_feces <- 0
  for (duration in random_duration_rounded) {
     prod_infected_pop_infected_feces<-prod_infected_pop_infected_feces+ duration*rate  #[g feces / year]
}

# Yearly production of non-infected feces by the infected population :
prod_infected_pop_non_infected_feces <- 0
for (duration in random_duration_rounded) {
  prod_infected_pop_non_infected_feces<-prod_infected_pop_non_infected_feces+ (365-duration)*rate  #[g feces / year]
}

#  Total yealy quantity of uninfected feces :
   tot_non_infected_feces <- prod_infected_pop_non_infected_feces+prod_no_infected_pop


   
# Density of microorganisms in the feces sample [OPT1] : 
  #Generate a fixed excretion density for each person

  # Set parameters for Density [log10 nb org / g feces]
    min_density <- 4   # Minimum value
    max_density <- 10   # Maximum value
    mode_density <- 6   # Mode (most probable value)
    n <- X*incidence   # Number of random numbers to generate (= infected population)

  # Generate random numbers from a triangular distribution : 
    # attribute a density of excretion [ nb org / g feces] to each infected person
    random_density_OPT1 <- ceiling(10^(rtriangle(n, min_density, max_density, mode_density)))
      # here we take 10^() to convert from log10(nb. of org.) to nb. of org.
      # also we rounded the value because we can't have half of an organism
  
   
    
  #Total number of microorganism in infected feces
   #Compute the concentration of microb. in the yearly production of infected feces : 
   nb_microb_in_infected_feces_OPT1 <- 0
   for (i in 1:length(random_duration_rounded)) {
     nb_microb_in_infected_feces_OPT1<-nb_microb_in_infected_feces_OPT1+ random_duration_rounded[i]*rate*random_density_OPT1[i]
    }


#(there's no OPT2)

   
# Density of microorganisms [OPT 3]: BEST ONE(?)
  # Generate an excretion density for each excretion of each person


  #Compute the concentration of microb. in the yearly production of infected feces : 
   nb_microb_in_infected_feces_OPT3 <- 0 #[nb of org in the yearly production of infected feces]

   for (duration in random_duration_rounded) { #go through each person
 
     # Generate random numbers from a triangular distribution : 
     # attribute a density of excretion [ nb org / g feces] to each infected excretion 
  
     # Set parameters for Density [log10 nb org / g feces]
        min_density <- 4   # Minimum value
        max_density <- 10   # Maximum value
        mode_density <- 6   # Mode (most probable value)
        n <-duration  # Number of random numbers to generate (= each time there's an infected excretion)
    
      random_density_OPT3 <- (10^(rtriangle(n, min_density, max_density, mode_density))) #[nb of org/g of fece]
        # here we take 10^() to convert from log10(nb. of org.) to nb. of org.
        # also we rounded the value because we can't have half of an organism
      
     for(density in random_density_OPT3){ #go through each excretion 
        nb_microb_in_infected_feces_OPT3<- nb_microb_in_infected_feces_OPT3+ rate*density 
  }
  }

# Let's now compare the two results :
print("OPT 1: ")
print(nb_microb_in_infected_feces_OPT1)
print("OPT 3: ")
print(nb_microb_in_infected_feces_OPT3)

pourcentage_difference <- ((nb_microb_in_infected_feces_OPT3 - nb_microb_in_infected_feces_OPT1) / nb_microb_in_infected_feces_OPT1) * 100
print(pourcentage_difference)


# Now let's compute the average density of the pathogen in mixed feces from a population :
   average_pathogen_density_OPT1<-nb_microb_in_infected_feces_OPT1/(tot_non_infected_feces+prod_infected_pop_infected_feces)
   average_pathogen_density_OPT3<- nb_microb_in_infected_feces_OPT3/(tot_non_infected_feces+prod_infected_pop_infected_feces)
  #[nb org excreted / year ÷ g feces/ year] => [nb. org.  / g feces in mixed population]
   
   # Let's now compare the two results :
   print("OPT 1: pathogen density in mixed feces sample [nb ogr / g sample] ")
   print(average_pathogen_density_OPT1)
   print("OPT 3: pathogen density in mixed feces sample [nb ogr / g sample]")
   print(average_pathogen_density_OPT3)
   
   print("OPT 1: pathogen density in infected feces  [nb ogr / g sample] ")
   print(nb_microb_in_infected_feces_OPT1/(prod_infected_pop_infected_feces))
   print("OPT 3: pathogen density in infected feces sample [nb ogr / g sample]")
   print(nb_microb_in_infected_feces_OPT3/(prod_infected_pop_infected_feces))
   
