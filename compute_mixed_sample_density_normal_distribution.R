#clear all variables :
rm(list = ls())

#####################HAZARD IDENTIFICATION##################################################

# create dataframe of concentration => made by the data of SCH�NNING ET AL, 2006 [ Ln of organisms / g of feces] 
# Those are numbers following a normal distribution

# Data for duration of excretion [days], ln
duration_excretion <- data.frame(
  salm = c(3.6, 0.2),
  rota = c(1.6, 1.23),
  crypto = c(2.0,0.85),
  giardia = c(4.5, 0.7),
  EHEC = c(2.1, 0.25)
)

# Data for Excretion density [nb of org/g of feces], ln
excretion_density <- data.frame(
  salm = c(13.8, 2.3),
  rota = c(20.7,2.3),
  crypto = c(17.3, 0.6),
  giardia = c(15.0,1.7),
  EHEC = c(5.8, 1.2)
)

# Data for incidences [ ]
inc_df <- data.frame(
  salm = c(500/100000,100/100000),
  rota = c(1200/100000,200/100000),
  crypto = c(200/100000, 25/100000),
  giardia = c(1100/100000,100/100000),
  EHEC = c(30/100000, 5/100000)
)

# Other parameters
X <- 500000             # fictive population size    
rate <- 250             # [g/cap/day] #Study of TO-FIGUERAS ET AL.(2000), BALCELLS GORINA (1989), SCHOUW ET AL.(2002)

pathogens <- c("salm", "rota", "crypto", "giardia", "EHEC")


# Initialize lists to store results for each pathogen
prod_infected_pop_infected_feces_list <- list()
prod_infected_pop_non_infected_feces_list <- list()


#####################Production of feces, infected and non infected############################

# Set parameters for Duration [days]
mean_duration <- duration_excretion[1,]   # Minimum value
sd_duration <- duration_excretion[2,] # Mode (most probable value)

n <- inc_df[1,]*X # this is a generalization: we took the mean without taking into consideration the standard deviation. this will be changed towards a distribution in a further step

prod_no_infected_pop <- list()
durations <- list()
for (pathogen in pathogens) {
  # Production of solid excreta of the non infected population => will be represented as well as a distribution on a further step
  prod_no_infected_pop[[pathogen]] <- mean(rate*365*(1-rnorm(X, mean = inc_df[1, pathogen], sd = inc_df[2,pathogen]))*X) # [g/year] 
  
  # Calculates the durations
  durations[[pathogen]] <- (rnorm(n[1, pathogen], mean_duration[1, pathogen], sd_duration[1, pathogen])) # [days/year]
  durations[[pathogen]] = ceiling(exp(durations[[pathogen]]))
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


#Calculation of total non infected feaces
tot_non_infected_feaces <- list()

for (pathogen in pathogens) {
  tot_non_infected_feaces[[pathogen]] <- prod_infected_pop_non_infected_feces_list[[pathogen]]+prod_no_infected_pop[[pathogen]]
}


#####################Density of microorganisms in the feaces#####################################

# Calculation of the number of pathogens in the infected feaces
nb_microb_in_infected_feces_list <- list()

mean_density <- excretion_density[1,]   # Mean
sd_density <- excretion_density[2,]   # Standard deviation

for (pathogen in names(durations)) {
  nb_microb_in_infected_feces <- 0
  for (duration in durations[[pathogen]]) {
    n <- duration
    random_density <- (exp(rnorm(n, mean_density[1, pathogen], sd_density[1, pathogen])))
    for (density in random_density) {
      nb_microb_in_infected_feces <- nb_microb_in_infected_feces + rate * density
    }
  }
  nb_microb_in_infected_feces_list[[pathogen]] <- nb_microb_in_infected_feces 
} # [nb of pathogens in the yearly production of infected feces]


# Find the average pathogen density of each pathogen
average_pathogen_density_list <- list()

for (pathogen in names(nb_microb_in_infected_feces_list)) {
  average_pathogen_density_list[[pathogen]] <- 
    nb_microb_in_infected_feces_list[[pathogen]] / (tot_non_infected_feaces[[pathogen]] + prod_infected_pop_infected_feces_list[[pathogen]])
} # [nb org excreted / year ÷ g feces/ year] => [nb. org.  / g feces in mixed population]



#####################EXPOSURE PATHWAY#####################################

# Transfer rates as indicated in Brooks: 
## The parameters here are an approximation of their distribution. 
## In the final code, this parameters will be in form of probability distribution
t_hm <- 0.36 #[]
t_fh <- 0.43 # []
t_gf <- 0.27 # []
t_exg <- (0.16+0.28)/2 # [mg/cm^2]
a_gex <- (0.13+0.25)/2*420 # [cm^2]

ingestion_g <- t_exg*a_gex*t_fh*t_hm*t_gf

# compute the number of pathogens ingested
n_pathogens_ingested <- list()

for (pathogen in names(average_pathogen_density_list)) {
  n_pathogens_ingested[[pathogen]] = 
    ingestion_g * average_pathogen_density_list[[pathogen]]/1000
} #[nb org excreted / year ÷ g feces/ year] => [nb. org.  / mg feces in mixed population]

#####################DOSE-RESPONSE###########################################################

pathogens <- c("salm", "rota", "crypto", "giardia", "EHEC")

# Exponential model
response_exponential_crypto <- 1-exp(-0.0042*n_pathogens_ingested[["crypto"]]) 
response_exponential_giardia <- 1-exp(-0.0198*n_pathogens_ingested[["giardia"]]) 

# Beta poisson
response_beta_poisson_rota = 1-(1+n_pathogens_ingested[["rota"]]/(0.42))^(-0.26) #Dose needs to be adjusted to the INGESTED DOSE
response_beta_poisson_salmonellanontyphoide <- 1-(1+n_pathogens_ingested[["salm"]]*(2^(1/(0.21))-1)/49.8)^(-0.21)
response_beta_poisson_EHEC <- 1-(1+n_pathogens_ingested[["EHEC"]]/(39.71))^(-0.373)


   

  
  
  

  
