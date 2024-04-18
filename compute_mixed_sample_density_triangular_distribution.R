#clear all variables :
  rm(list = ls())

#library for traingular distribution :
  library(triangle)
  library(EnvStats)

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
    adeno = c(8,10,12),
    rota = c(8,10,12),
    campylo = c(4,6,10),
    salm = c(6,7.5,9),
    crypto = c(6,7,9)
  )

# Data for incidences (same soruce as trinangular distribution) [ ]
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

pathogens <- c("adeno", "rota", "campylo", "salm", "crypto")


# Initialize lists to store results for each pathogen
prod_infected_pop_infected_feces_list <- list()
prod_infected_pop_non_infected_feces_list <- list()


#####################Production of feces, infected and non infected############################

# Yearly production of non-infected feces by the non infected population :
prod_no_infected_pop <- rate*365*(1-inc_df)*X # [g/year] 

# Yearly production of infected feces by the infected population :
# Set parameters for Duration [days]
min_duration <- duration_excretion[1,]   # Minimum value
max_duration <- duration_excretion[3,]  # Maximum value
mode_duration <- duration_excretion[2,]  # Mode (most probable value)
n <- inc_df[1,]*X  # Number of random numbers to generate (= infected population)


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
    random_density <- (10^(rtriangle(n, min_density[1, pathogen], max_density[1, pathogen], mode_density[1, pathogen])))
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
# The multiplication by 10^3 is to go from g to mg 





#####################EXPOSURE PATHWAY#########################################################

#Transfer rates as indicated in Brooks
t_hm <- 0.36 #[]
t_fh <- 0.43
t_gf <- 0.27 #to be changed
#a_gf <- (0.1+0.17)/2*420
t_exg <- (0.16+0.28)/2 #[mg/cm^2]
a_gex <- (0.13+0.25)/2*420 #[cm^2]

#g ingested
ingestion_g <- t_exg*a_gex*t_fh*t_hm*t_gf

# N pathogens ingested
n_pathogens_ingested <- list()

# Loop through each pathogen
for (pathogen in names(average_pathogen_density_list)) {
  # Perform element-wise division and store the result in the result_list
  n_pathogens_ingested[[pathogen]] = ingestion_g * average_pathogen_density_list[[pathogen]]/1000
  print('hi')
} #[nb org excreted / year ÷ g feces/ year] => [nb. org.  / mg feces in mixed population]
# division by 1000 is to go from g to mg in order to have consistent units. 

#something dont work here. it shoulld diminish by 10^(-5), no? whats wrong?

#####################DOSE-RESPONSE###########################################################

# Data for dose response model
param_beta_poisson <- data.frame(
  rota = c(0.26, 0.42), # c(alpha, beta)
  campylo1996 = c(0.145, 7.59),
  campylo2005 = c(0.024, 0.011)
)

param_exponential <- data.frame(
  adeno = 4.172*10^(-1),
  crypto = 0.0042
)



# Exponential model
response_exponential_adeno = 1-exp(-4.172*10^(-1)*n_pathogens_ingested[["adeno"]]) #Dose needs to be adjusted to the INGESTED DOSE
response_exponential_crypto = 1-exp(-0.0042*n_pathogens_ingested[["crypto"]]) 

# Beta poisson
response_beta_poisson_rota = 1-(1+n_pathogens_ingested[["rota"]]/(0.26))^(-0.42) #Dose needs to be adjusted to the INGESTED DOSE
response_beta_poisson_campylo1996 = 1-(1+n_pathogens_ingested[["campylo"]]/(0.145))^(-7.59) #Campylo jejuni
response_beta_poisson_campylo2005 = 1-(1+n_pathogens_ingested[["campylo"]]/(0.024))^(-0.011)
response_beta_poisson_salmonellanontyphoide = 1-(1+n_pathogens_ingested[["salm"]]*(2^(1/(0.21))-1)/49.8)^(-0.21)
###########!!!!!!!!!!!!!!!Attention, samonella probablement not real route, check!!!!!!!!!!!
#####################RISK CHARACTERIZATION##################################################


  
  
  

  
