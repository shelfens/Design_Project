#clear all variables :
  rm(list = ls())

#library for traingular distribution :
  library(triangle)
  library(EnvStats)

# The data given by Sch�nning follows a normal distribution.
  
# We therefore create samples of 10 datapoints following the given normal distribution per pathogen
 # conc_normal_dist <- data.frame(
  #  rota = c(20.7,2.3,9.8*10^8),
   # ecoli = c(5.8,1.2, 3.3*10^2),
    #salmonella = c(13.8,2.3,9.9*10^5),
    #crypto =  c(17.3,0.6,3.3*10^7),
    #giardia = c(15.0,1.7,3.3*10^6)
#  )

#####################HAZARD IDENTIFICATION##################################################
#####################Definition of parameters###############################################
  
# create dataframe of concentration => made by the data of PETTERSON ET AL, 2016 [ Log10 of organisms / g of feces] 
# Those are numbers following a triangular distribution => needs to be adapted
  
# Data for duration of excretion [days]
  duration_excretion <- data.frame(
    adeno = c(3,7,12),
    rota = c(3,7,12), 
    campylo = c(15,34,42),
    salm = c(10,15,50),
    crypto = c(5,10,30)
  )
  
# Data for Excretion density [nb of org/g of feces]
  excretion_density <- data.frame(
    adeno = c(10 ^ c(8,10,12)),
    rota = c(10^c(8,10,12)),
    campylo = c(10 ^ c(4,6,10)),
    salm = c(10 ^ c(6,7.5,9)),
    crypto = c(10 ^ c(6,7,9))
  )

# Data for incidences (same soruce as trinangular distribution) [%]
  inc_df <- data.frame(
    adeno = 0.88/100,
    rota = 1.28/100,
    campylo = 2.4/100,
    salm = 0.88/100,
    crypto = 1.28/100
  )

# Other parameters
X <- 500000             # fictive population size    
rate <- 250             # [g/cap/day] 

pathogens <- c("rota", "adeno", "campylo", "salm", "crypto")


# Initialize lists to store results for each pathogen
prod_infected_pop_infected_feces_list <- list()
prod_infected_pop_non_infected_feces_list <- list()

# Loop through each pathogen
for (pathogen in pathogens) {
  # Set parameters for Duration [days]
  min_duration <- duration_excretion[1, pathogen]   # Minimum value
  max_duration <- duration_excretion[3, pathogen]  # Maximum value
  mode_duration <- duration_excretion[2, pathogen]  # Mode (most probable value)
  n <- inc_df[1, pathogen] * X  # Number of random numbers to generate (= infected population)
  
  # Generate random numbers from a triangular distribution:
  # attribute a duration of excretion [days] to each infected person
  # We need to round to nearest whole number the values of duration
  # (to have entire days)
  random_duration <- rtriangle(n, min_duration, max_duration, mode_duration)
  random_duration_rounded <- ceiling(random_duration)
  
  # Compute the yearly production of infected feces by the infected population:
  prod_infected_pop_infected_feces <- sum(random_duration_rounded * rate)  # [g feces / year]
  
  # Compute the yearly production of non-infected feces by the infected population:
  prod_infected_pop_non_infected_feces <- sum((365 - random_duration_rounded) * rate)  # [g feces / year]
  
  # Store results for this pathogen
  prod_infected_pop_infected_feces_list[[pathogen]] <- prod_infected_pop_infected_feces
  prod_infected_pop_non_infected_feces_list[[pathogen]] <- prod_infected_pop_non_infected_feces
}




#####################Production of feces, infected and non infected############################

# Yearly production of non-infected feces by the non infected population :
prod_no_infected_pop <- 250*365*(1-inc_df)*X # [g/year] 

# Yearly production of infected feces by the infected population :
# Set parameters for Duration [days]
min_duration <- duration_excretion[1,]   # Minimum value
max_duration <- duration_excretion[3,]  # Maximum value
mode_duration <- duration_excretion[2,]  # Mode (most probable value)
n <- inc_df[1,]*X  # Number of random numbers to generate (= infected population)

# Generate random numbers from a triangular distribution : 
# attribute a duration of excretion [days] to each infected person
# We need to round to nearest whole number the values of duration
# (to have entires days)
durations <- list(
  rota = ceiling(rtriangle(n[1,'rota'], min_duration[1, 'rota'], max_duration[1, 'rota'], mode_duration[1,'rota'])),
  adeno = ceiling(rtriangle(n[1,'adeno'], min_duration[1, 'adeno'], max_duration[1, 'adeno'], mode_duration[1,'adeno'])),
  campylo = ceiling(rtriangle(n[1,'campylo'], min_duration[1, 'campylo'], max_duration[1, 'campylo'], mode_duration[1,'campylo'])),
  salm = ceiling(rtriangle(n[1,'salm'], min_duration[1, 'salm'], max_duration[1, 'salm'], mode_duration[1,'salm'])),
  crypto = ceiling(rtriangle(n[1,'crypto'], min_duration[1, 'crypto'], max_duration[1, 'crypto'], mode_duration[1,'crypto']))
)

prod_infected_pop_infected_feces_list <- list()
prod_infected_pop_non_infected_feces_list <- list()

for (pathogen in names(durations)) {
  prod_infected_pop_infected_feces <- 0
  prod_infected_pop_non_infected_feces <- 0
  for (duration in durations[[pathogen]]) {
    prod_infected_pop_infected_feces <- prod_infected_pop_infected_feces + duration * rate
    prod_infected_pop_non_infected_feces<-prod_infected_pop_non_infected_feces+ (365-duration)*rate
  }
  prod_infected_pop_infected_feces_list[[pathogen]] <- prod_infected_pop_infected_feces
  prod_infected_pop_non_infected_feces_list[[pathogen]] <- prod_infected_pop_non_infected_feces
}


#  Total yealy quantity of uninfected feces :
tot_non_infected_feces <- prod_infected_pop_non_infected_feces+prod_no_infected_pop


#####################Density of microorganisms in the feaces#####################################

# Generate an excretion density for each excretion of each person

nb_microb_in_infected_feces_list <- list()
   
min_density <- excretion_density[1,]   # Minimum value
max_density <- excretion_density[3,]  # Maximum value
mode_density <- excretion_density[2,]  # Mode (most probable value)

#####Not that sure, that it is working the right way. To check. ####
for (pathogen in names(durations)) { #for each pathogen
  nb_microb_in_infected_feces <- 0
  for (duration in durations[[pathogen]]) {  # Change "campylo" to [pathogen]
    n <- duration
    random_density <- rtriangle(n, min_density[1, pathogen], max_density[1, pathogen], mode_density[1, pathogen])
    for (density in random_density) {
      nb_microb_in_infected_feces <- nb_microb_in_infected_feces + rate * density
    }
  }
  nb_microb_in_infected_feces_list[[pathogen]] <- nb_microb_in_infected_feces
  print(nb_microb_in_infected_feces)#[nb of org in the yearly production of infected feces]
}



# Find the average pathogen density of each pathogen
average_pathogen_density_list <- list()

# Loop through each pathogen
for (pathogen in names(nb_microb_in_infected_feces_list)) {
  # Perform element-wise division and store the result in the result_list
  average_pathogen_density_list[[pathogen]] <- nb_microb_in_infected_feces_list[[pathogen]] / (tot_non_infected_feces[[pathogen]] + prod_infected_pop_infected_feces_list[[pathogen]])
} #[nb org excreted / year ÷ g feces/ year] => [nb. org.  / g feces in mixed population]





#####################EXPOSURE PATHWAY#########################################################
ingestion <- 0.1 # Assuming the 0.1g ingested as done in the base study (Brooks)

#Transfer rates as indicated in Brooks
fomite_to_hand <- 0.43
hand_to_mouth <-0.36

#g ingested
ingestion_g <- 0.1*fomite_to_hand*hand_to_mouth

# N pathogens ingested
n_pathogens_ingested <- list()

# Loop through each pathogen
for (pathogen in names(average_pathogen_density_list)) {
  # Perform element-wise division and store the result in the result_list
  n_pathogens_ingested[[pathogen]] <- ingestion_g * average_pathogen_density_list[[pathogen]]
} #[nb org excreted / year ÷ g feces/ year] => [nb. org.  / g feces in mixed population]




#####################DOSE-RESPONSE###########################################################

# Data for dose response model
param_beta_poisson <- data.frame(
  rota = c(0.26, 0.42) # c(alpha, beta)
)

param_exponential <- data.frame(
  adeno = 4.172*10^(-1),
  crypto = 0.0042
)


response_exponential_adeno = 1-exp(-4.172*10^(-1)*n_pathogens_ingested[["adeno"]]) #Dose needs to be adjusted to the INGESTED DOSE

response_beta_poisson = 1-(1+dose/(param_beta_poisson[2,])^(-param_beta_poisson[1,])) #Dose needs to be adjusted to the INGESTED DOSE



#####################RISK CHARACTERIZATION##################################################





################################to be continued##############################################
   # Let's now compare the two results :
   print("pathogen density in mixed feces sample [nb ogr / g sample]")
   print(average_pathogen_density_list)
   
   
   
   
   
   
  
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
           

            
      
  
   

  
  
  

  
