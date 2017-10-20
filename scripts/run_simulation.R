library(reshape2)
library(ggplot2)

wd <- "~/Global Burden/vaxedemic/" ## working directory ##DH - deleted "Documents/" - R is weird!
source(paste0(wd, "R/simulation.R"))
source(paste0(wd, "R/setup.R"))

## LIFE HISTORY PARAMETER INPUTS
## R_0, recovery time and latent period
life_history_params <- list(R0=1.1, TR=2.6, LP = 1.5)

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
                         country_specific_contact = FALSE,
                         seasonal = FALSE,
                         rng_seed = 0)

## run simulation for tmax days
tmax <- 100
## tdiv timesteps per day
tdiv <- 24 ##DH
## allocate and distribute vaccine every vac_alloc_period time divisions
## i.e. in this example, every 7 days
vax_alloc_period <- 24 * 7

# set random number generation seed
if(!is.null(simulation_flags[["rng_seed"]])) {
  set.seed(simulation_flags[["rng_seed"]])
}

if(simulation_flags[["real_data"]]) {
  ## get number of countries and ages from files
  demography_filename <- paste0(wd, "data/unified/demographic_data_intersect.csv")
  contact_filename <- paste0(wd, "data/unified/contact_data_intersect.csv")

  travel_filename <- paste0(wd, "data/unified/flight_data_intersect.csv")
  latitude_filename <- paste0(wd, "data/unified/latitudes_intersect.csv")
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
if(simulation_flags[["real_data"]]) {
  seedCountries <- "China"
} else {
  seedCountries <- 1
}

## number of exposed individuals in each seeded country
seedSizes <- c(20)
## which age group(s) to seed
seedAges <- 3
## which risk group(s) to seed
seedRiskGroups <- 1

if(simulation_flags[["real_data"]]) {
  ## construct demography matrix
  tmp <- setup_populations_real_data(demography_filename,
                            risk_propns, risk_factors)

  ## construct travel matrix
  K <- setup_travel_real_data(travel_filename, tmp$pop_size, travel_params)
  ## construct latitude vector
  latitudes <- read_latitude_data(latitude_filename)

} else {
  
  ## construct demography matrix
  tmp <- setup_populations(popn_size,n_countries,age_propns, 
                           risk_propns, risk_factors)
  ## construct travel matrix
  K <- matrix(1,n_countries,n_countries)+999*diag(n_countries) #Travel coupling - assumed independent of age (but can be changed)
  ## construct latitude vector
  latitudes <- matrix(seq_len(n_countries), n_countries, 1)
}

X <- tmp$X
labels <- tmp$labels

#construct vector of number of exposed individuals in each location, age, risk group 
## for now, can only seed in one location, age, risk group. Vectorise later
seed_vec <- double(length(X))
seed_vec[(which(labels$Location == seedCountries & labels$Age == seedAges &
              labels$RiskGroup == seedRiskGroups))[1]] <- seedSizes

## Generate a contact matrix with dimensions (n_ages*n_riskgroups) * (n_ages*n_riskgroups). ie. get age specific,
## then enumerate out by risk group. If we have country specific contact rates, we get a list
## of these matrices of length n_countries
if(simulation_flags[["real_data"]]) {
  C1 <- read_contact_data(contact_filename)
} else {
  ## Contact rates
  contactRates <- c(6.92,.25,.77,.45,.19,3.51,.57,.2,.42,.38,1.4,.17,.36,.44,1.03,1.83)
  contactDur <- c(3.88,.28,1.04,.49,.53,2.51,.75,.5,1.31,.8,1.14,.47,1,.85,.88,1.73)
  # for now, make contact matrices same for all countries
  C1 <- generate_contact_matrix(contactRates, contactDur,n_ages, simulation_flags[["ageMixing"]])
  if(simulation_flags[["country_specific_contact"]]) {
    C1 <- rep(list(C1), n_countries)
  }
}

## Generate risk factor modifier. ie. modifier for each age/risk group pair, same dimensions as C2
## The risk factor modifier modifies the susceptibility of age/risk groups
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
  ## vaccine production function defined here
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

## make the vaccine production function using the above closure
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
  ## vaccine allocation function defined here
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

# make the vaccine alloation function using the above closure
vax_allocation_func <- vaccine_allocation_closure(K, vax_allocation_params, labels)

## gather simulation parameters
sim_params <- list(n_countries=n_countries,
                   n_ages=n_ages,
                   n_riskgroups=n_riskgroups,
                   seed_vec = seed_vec)

## run simulation
res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                      X, C3, K, latitudes, cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period)

## plot stuff 
plot_labels <- expand.grid("Time"=seq(0,tmax,by=1/tdiv),"Location"=1:n_countries,"Age"=1:n_ages,"RiskGroup"=1:n_riskgroups)

tend <- tmax*tdiv+1
popTotal=sum(labels[,1])
globalAttack <- sum(res$R[,tend] + res$RV[,tend])/popTotal

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
