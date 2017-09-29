library(reshape2)
library(ggplot2)
source("~/Documents/vaxedemic/R/simulation.R")
source("~/Documents/vaxedemic/R/setup.R")

## LIFE HISTORY PARAMETER INPUTS
## R_0, recovery time and latent period
life_history_params <- list(R0=1.8, TR=2.6, LP = 1.5)

## travel parameters: scaling of off-diagonals
travel_params <- list(epsilon = 1e-3)

## vaccine efficacy and initial vaccinated proportion
# this example roughly brings effective R to 1.2
vax_params <- list(efficacy = 1 - 1.2/1.8, propn_vax0 = 0)

## example parameters for vaccine production (see cum_vax_pool_func_closure)
vax_production_params <- list(detection_delay = 0, production_delay = 0, 
                              production_rate = 0, max_vax = 5e9)

## example parameters for vaccine allocation
vax_allocation_params <- list(example = 1)

## SIMULATION OPTIONS
simulation_flags <- list(ageMixing=TRUE,
                         riskGroups=TRUE,
                         normaliseTravel=TRUE,
                         spatialCoupling=TRUE,
                         real_data = TRUE,
                         country_specific_contact = TRUE,
                         seed = 1)
tmax <- 100
tdiv <- 24
## allocate and distribute vaccine every vac_alloc_period time divisions
## i.e. in this example, every 7 days
vax_alloc_period <- 24 * 7 

if(!is.null(simulation_flags[["seed"]])) {
  set.seed(simulation_flags[["seed"]])
}

if(simulation_flags[["real_data"]]) {
  demography_filename <- "~/Documents/vaxedemic/data/unified/demographic_data_intersect.csv"
  contact_filename <- "~/Documents/vaxedemic/data/unified/contact_data_intersect.csv"
  travel_filename <- "~/Documents/vaxedemic/data/unified/flight_data_intersect.csv"
  tmp <- read.csv(demography_filename, sep = ",")
  n_countries <- nrow(tmp)
  n_ages <- ncol(tmp) - 2
} else {
  ## SETUP FAKE COUNTRY DATA
  popn_size <- 100000
  
  
  ## Setup age propns
  n_ages <- 4
  age_propns <- rep(1/n_ages, n_ages)
  age_propns <- c(5,14,45,16)/80
  n_countries <- 10
}

## Setup risk groups
n_riskgroups <- 2

risk_propns <- rep(1/n_riskgroups,n_riskgroups) ## Assume risk groups are uniformly distributed
risk_propns <- matrix(rep(risk_propns,each=n_ages),ncol=n_riskgroups) ## Assume that proportion of ages in each risk group are the same for all ages
risk_factors <- rep(1, n_riskgroups) ## Assume that each risk group has same modifier

## Enumerate out risk factors for each age group
age_specific_riskgroup_factors <- matrix(rep(risk_factors,each=n_ages),
                                         ncol=n_riskgroups)

## Seeding setting
seedCountries <- c(1)
seedSizes <- c(10)
seedAges <- 3

## Contact rates
contactRates <- c(6.92,.25,.77,.45,.19,3.51,.57,.2,.42,.38,1.4,.17,.36,.44,1.03,1.83)
contactDur <- c(3.88,.28,1.04,.49,.53,2.51,.75,.5,1.31,.8,1.14,.47,1,.85,.88,1.73)



## to do: if(simulation_flags[["real_data"]]), read in real data 


if(simulation_flags[["real_data"]]) {
  ## demography
  tmp <- setup_populations_real_data(demography_filename,
                            risk_propns, risk_factors,
                            n_riskgroups)
  ## Travel coupling
  K <- setup_travel_real_data(travel_filename, tmp$pop_size, travel_params)
} else {c
  
  ## demography
  tmp <- setup_populations(popn_size,n_countries,age_propns, n_ages,
                           risk_propns, risk_factors,
                           n_riskgroups)
  ## Travel coupling
  K <- matrix(1,n_countries,n_countries)+999*diag(n_countries) #Travel coupling - assumed independent of age (but can be changed)
}

X <- tmp$X
labels <- tmp$labels

## Generate a contact matrix with dimensions (n_ages*n_riskgroups) * (n_ages*n_riskgroups). ie. get age specific,
## then enumerate out by risk group. If we had country specific contact rates, we get a list
## of these matrices of length n_countries
if(simulation_flags[["real_data"]]) {
  C1 <- read_contact_data(contact_filename)
} else {
  C1 <- generate_contact_matrix(contactRates, contactDur,n_ages, simulation_flags[["ageMixing"]])
  if(simulation_flags[["country_specific_contact"]]) {
    C1 <- rep(list(C1), n_countries)
  }
}


## Generate risk factor modifier. ie. modifier for each age/risk group pair, same dimensions as C2
risk <- c(t(age_specific_riskgroup_factors))
risk_matrix <- t(kronecker(risk,matrix(1,1,n_riskgroups*n_ages)))
if(is.list(C1)) { # for country specific contact rates
  C2 <- lapply(C1, function(x) kronecker(x, matrix(1,n_riskgroups,n_riskgroups)))
  C3 <- lapply(C2, function(x) x*risk_matrix)
} else {
  C2 <- kronecker(C1, matrix(1,n_riskgroups,n_riskgroups))
  C3 <- C2*risk_matrix
}


## vaccination production curve

# a closure to make a function which returns a scalar:
# the number of vaccines which have ever existed at time t
# (may involve integrating the production rate curve if you want,
# otherwise provide directly)
# current example:
# no vaccine produced until time vax_production_params[["detection_delay"]] + 
# vax_production_params[["production_delay"]], then
# constant production rate until max number of doses ever made reached,
# then no production

cum_vax_pool_func_closure <- function(vax_production_params) {
  function(t) {
    t_since_production <- t - (vax_production_params[["detection_delay"]] + 
      vax_production_params[["production_delay"]])
    if(t_since_production < 0) {
      0
    } else {
      min(vax_production_params[["max_vax"]],
          t_since_production * vax_production_params[["production_rate"]])
    }
  }
}

cum_vax_pool_func <- cum_vax_pool_func_closure(vax_production_params)

## vaccine allocation strategy

# a closure to make a function which returns a vector of length
# n_countries * n_ages * n_risk_groups:
# the number of vaccines to be allocated to each location, age, risk group
# vaccine allocation is a function of the current state of the epidemic, the
# rate of change of the epidemic (i.e. incidence etc.),
# current vaccine pool size, the travel matrix and other parameters

# probably can make independent of the vaccinated classes...

# if non-integer outputted, main simulation will round the numbers down
# this may cause problems because of the number of classes, might be that
# no vaccines are ever allocated -- we'll see

# current example:
# allocate vaccines according to absolute incidence in each location/age/risk group
# note unrealism because we would obviously want to allocate vaccines to everyone 
# in the same country as the epidemic is occurring regardless of the age/risk group
# in which the current infectious people are

# to improve this, I've included the labels as an input for ease of coding, so that
# we can easily figure out the elements in each vector corresponding to
# e.g. the people in the same country as people currently infected
vaccine_allocation_closure <- function(travel_matrix, vax_allocation_params, labels) {

  function(S, E, I, R, SV, EV, IV, RV, vax_pool) {
    if(any(E > 0)) {
      return(E / sum(E) * vax_pool) # incidence proportional to E
    } else if(any(I > 0)){
      return(I / sum(I) * vax_pool) # if no exposed unvaccinated left, allocate proportional to prevalence
    } else if(any(S > 0)){
      return(S / sum(S) * vax_pool) # if no infectious unvaccinated left, allocate proportional to susceptibles
    } else {
      return(S * 0) # allocate nothing
    }
  }
}

vax_allocation_func <- vaccine_allocation_closure(K, vax_allocation_params, labels)

## Normalise 

sim_params <- list(n_countries=n_countries,
                   n_ages=n_ages,
                   n_riskgroups=n_riskgroups,
                   seedCs=seedCountries,
                   seedNs=seedSizes,
                   seedAges=seedAges)

res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                      X, C3, K, cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period)

plot_labels <- expand.grid("Time"=seq(0,tmax,by=1/tdiv),"Location"=1:n_countries,"Age"=1:n_ages,"RiskGroup"=1:n_riskgroups)

I <- cbind(labels[,c("Location","Age","RiskGroup")], res$I + res$IV)
I <- melt(I, id.vars=c("Location","Age","RiskGroup"))
I$Age <- as.factor(I$Age)
I$RiskGroup <- as.factor(I$RiskGroup)
I$variable <- as.numeric(I$variable)
times <- seq(0,tmax,by=1/tdiv)
I$variable <- times[I$variable]
I_aggregated <- aggregate(I[,"value"], I[,c("variable","Location","Age")], FUN=sum)

N <- aggregate(data=labels, X~Location + Age,FUN=sum)
I_aggregated <- merge(I_aggregated,N,id.vars=c("Location","Age"))

p1 <- ggplot(I_aggregated,aes(x=variable,y=x/X,col=Age)) +
    geom_line() +
    facet_wrap(~Location) +
    theme_bw()

# p2 <- ggplot(I, aes(x=variable,y=value,col=RiskGroup)) + geom_line() + facet_grid(Age~Location) + theme_bw()

# to do: make plot of vaccines allocated (SV + EV + IV + RV) over time by country
# grid_plot <- cowplot::plot_grid(p1,p2,ncol=2,align="hv")
