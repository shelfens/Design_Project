#clear all variables :
  rm(list = ls())

#library for traingular distribution :
  library(triangle)
  library(EnvStats)

#####################HAZARD IDENTIFICATION##################################################

# create dataframe of concentration => made by the data of PETTERSON ET AL, 2016 [ Log10 of organisms / g of feces] 
# Those are numbers following a triangular distribution
  
# Data for duration of excretion [days]
  duration_excretion <- data.frame(
    adeno = c(3,7,12),
    rota = c(3,7,12), 
    campylo = c(15,34,42),
    salm = c(10,15,50),
    crypto = c(5,10,30)
  )
  
# Data for Excretion density [nb of org/g of feces], log10
  excretion_density <- data.frame(
    adeno = c(8,10,12),
    rota = c(8,10,12),
    campylo = c(4,6,10),
    salm = c(6,7.5,9),
    crypto = c(6,7,9)
  )

# Data for incidences (same source as trinangular distribution) [ ]
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


#############Production of feces, infected and non infected##################

# Yearly production of non-infected feces by the non infected population :
prod_no_infected_pop <- rate*365*(1-inc_df)*X # [g/year] 

# Yearly production of infected feces by the infected population :
min_duration <- duration_excretion[1,]   # Minimum value
max_duration <- duration_excretion[3,]  # Maximum value
mode_duration <- duration_excretion[2,]  # Mode (most probable value)
n <- inc_df[1,]*X  # Number of random numbers to generate (= infected population)

durations <- list()
for (pathogen in pathogens) {
  # Calculates the durations
  durations[[pathogen]] <- (rtriangle(n[1, pathogen], min_duration[1, pathogen], max_duration[1, pathogen], mode_duration[1,pathogen])) # [days/year]
  durations[[pathogen]] = ceiling(durations[[pathogen]])
}


# Yearly production of infected and non-infected feces by the infected population :
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



################Density of microorganisms in the feaces#######################

# Compute the number of pathogens in the infected feaces
nb_microb_in_infected_feces_list <- list()
   
min_density <- excretion_density[1,]   # Minimum value
max_density <- excretion_density[3,]  # Maximum value
mode_density <- excretion_density[2,]  # Mode (most probable value)

for (pathogen in names(durations)) { 
  nb_microb_in_infected_feces <- 0
  for (duration in durations[[pathogen]]) {  
    n <- duration
    random_density <- (10^(rtriangle(n, min_density[1, pathogen], max_density[1, pathogen], mode_density[1, pathogen])))
    for (density in random_density) {
      nb_microb_in_infected_feces <- nb_microb_in_infected_feces + rate * density
    }
  }
  nb_microb_in_infected_feces_list[[pathogen]] <- nb_microb_in_infected_feces
} # [nb of org in the yearly production of infected feces]


# Find the average pathogen density of each pathogen
average_pathogen_density_list <- list()

for (pathogen in names(nb_microb_in_infected_feces_list)) {
  average_pathogen_density_list[[pathogen]] <- nb_microb_in_infected_feces_list[[pathogen]] / (tot_non_infected_feces[[pathogen]] + prod_infected_pop_infected_feces_list[[pathogen]])
} #[nb org excreted / year ÷ g feces/ year] => [nb. org.  / g feces in mixed population]


#####################EXPOSURE PATHWAY#########################################################

#Transfer rates as indicated in Brooks
t_hm <- 0.36 #[]
t_fh <- 0.43 #[]
t_gf <- 0.27 #[]
t_exg <- (0.16+0.28)/2 #[mg/cm^2]
a_gex <- (0.13+0.25)/2*420 #[cm^2]

#g ingested
ingestion_g <- t_exg*a_gex*t_fh*t_hm*t_gf

# N pathogens ingested
n_pathogens_ingested <- list()

for (pathogen in names(average_pathogen_density_list)) {
  n_pathogens_ingested[[pathogen]] = ingestion_g * average_pathogen_density_list[[pathogen]]/1000
} #[nb org excreted / year ÷ g feces/ year] => [nb. org.  / mg feces in mixed population]


#####################DOSE-RESPONSE###########################################################

# Exponential model
response_exponential_adeno = 1-exp(-4.172*10^(-1)*n_pathogens_ingested[["adeno"]]) #Dose needs to be adjusted to the INGESTED DOSE
response_exponential_crypto = 1-exp(-0.0042*n_pathogens_ingested[["crypto"]]) 

# Beta poisson
response_beta_poisson_rota = 1-(1+n_pathogens_ingested[["rota"]]/(0.42))^(-0.26) #Dose needs to be adjusted to the INGESTED DOSE
response_beta_poisson_campylo1996 = 1-(1+n_pathogens_ingested[["campylo"]]/(7.59))^(-0.145) #Campylo jejuni
response_beta_poisson_campylo2005 = 1-(1+n_pathogens_ingested[["campylo"]]/(0.011))^(-0.024)
  