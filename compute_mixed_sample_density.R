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


   
# Density of microorganisms :
  # Generate an excretion density for each excretion of each person


  #Compute the concentration of microb. in the yearly production of infected feces : 
   nb_microb_in_infected_feces <- 0 #[nb of org in the yearly production of infected feces]

   for (duration in random_duration_rounded) { #go through each person
 
     # Generate random numbers from a triangular distribution : 
     # attribute a density of excretion [ nb org / g feces] to each infected excretion 
  
     # Set parameters for Density [log10 nb org / g feces]
        min_density <- 4   # Minimum value
        max_density <- 10   # Maximum value
        mode_density <- 6   # Mode (most probable value)
        n <-duration  # Number of random numbers to generate (= each time there's an infected excretion)
    
      random_density <- (10^(rtriangle(n, min_density, max_density, mode_density))) #[nb of org/g of fece]
        # here we take 10^() to convert from log10(nb. of org.) to nb. of org.
        # also we rounded the value because we can't have half of an organism
      
     for(density in random_density){ #go through each excretion 
        nb_microb_in_infected_feces<- nb_microb_in_infected_feces+ rate*density 
  }
  }



# Now let's compute the average density of the pathogen in mixed feces from a population :
   average_pathogen_density<- nb_microb_in_infected_feces/(tot_non_infected_feces+prod_infected_pop_infected_feces)
  #[nb org excreted / year ÷ g feces/ year] => [nb. org.  / g feces in mixed population]
   
   # Let's now compare the two results :
   print("pathogen density in mixed feces sample [nb ogr / g sample]")
   print(average_pathogen_density)
   
   print("OPT 3: pathogen density in raw infected feces sample [nb ogr / g sample]")
   print(nb_microb_in_infected_feces/(prod_infected_pop_infected_feces))
   
  
    ################################################################
   
   
   #EXAMPLE FOR Salmonella ,study of Schonning et al. 2006: 
   
     X <- 500000             # fictive population size 
     rate <- 250             # [g/cap/day] 
     years <- 10             # nb of year to generate
     
    # Set parameters for Incidence 
     mean_incidence <- 500   # mean value
     sd_incidence <- 100     # standard deviation value
     years <-years  # Number of random numbers to generate (= one for each year)
     
     random_incidence <- rnorm(years, mean_incidence, sd_incidence)/1000 #[%]
       # here we divided per 1000 to have an incidence in [%] and not in per 100 000

    for (year in 1:years) {                # compute for each year with a different incidence
        
      incidence <-random_incidence[year]    # [%], select the incidence for the corresponding year
      
      # Yearly production of non-infected feces by the non infected population :
         prod_no_infected_pop <- 250*365*(1-incidence)*X  # [g/year]  
      
      # Yearly production of infected feces by the infected population :
         # Set parameters for Excretion Time [days]
         mean_duration <- 3.6   # mean value
         sd_duration <- 0.2     # standard deviation value
         median_duration <- 37  # median value
         n <- X*incidence  # Number of random numbers to generate (= infected population)
         
         random_duration <- rnorm(n, mean_duration, sd_duration) # TO MODIFY [days]
         random_duration_rounded = ceiling(random_duration)  
           # We need to round to nearest whole number the values of duration
           # (to have entires days)
        
           prod_infected_pop_infected_feces <- 0
           prod_infected_pop_non_infected_feces <- 0
           nb_microb_in_infected_feces <- 0 #[nb of org in the yearly production of infected feces]
           
           for (duration in random_duration_rounded) {  #go through each person
             #Compute the yearly production of infected and uninfected feces by the infected population : 
               prod_infected_pop_infected_feces<-prod_infected_pop_infected_feces+ duration*rate  #[g feces / year]
               prod_infected_pop_non_infected_feces<-prod_infected_pop_non_infected_feces+ (365-duration)*rate  #[g feces / year]
            
             #Compute the concentration of microb. in the yearly production of infected feces : 
               # Generate random numbers from a normal distribution : 
               # attribute a density of excretion [ nb org / g feces] to each infected excretion 
               
               # Set parameters for Density [log10 nb org / g feces]
               mean_density <- 13.8   # mean value
               sd_density <- 2.3   # standard deviation value
               median_duration <- 9.9*10^5   # median 
               n <-duration  # Number of random numbers to generate (= each time there's an infected excretion)
               
               random_density <- rnorm(n, mean_density, sd_density) # TO MODIFY 
               random_density_rounded = ceiling(random_density)
                 #[nb of org/g of fece]
                 # also we rounded the value because we can't have half of an organism
               
               for(density in random_density){ #go through each excretion 
                 nb_microb_in_infected_feces<- nb_microb_in_infected_feces+ rate*density 
               }
           }
           
                }
           

            
      
  
   

  
  
  

  
