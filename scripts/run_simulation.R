library(reshape2)
library(ggplot2)
library(Matrix)
library(data.table)

wd <- "~/Documents/vaxedemic/" 
devtools::load_all(wd)

## LIFE HISTORY PARAMETER INPUTS
## R_0, recovery time and latent period
life_history_params <- list(R0=1.8, TR=2.6, LP = 1.5, case_fatality_ratio = c(1e-3,1e-2))

## travel parameters: scaling of off-diagonals
travel_params <- list(epsilon = 1e-3)

## vaccine efficacy and initial vaccinated proportion
# this example roughly brings effective R to 1.2
vax_params <- list(efficacy = 1 - 1.2/1.8, propn_vax0 = 0)

## example parameters for vaccine production (see cum_vax_pool_func_closure)
vax_production_params <- list(detection_delay = 0, production_delay = 0, 
                              production_rate = 0, max_vax = 5e9)

## example parameters for vaccine allocation
vax_allocation_params <- list(priorities = NULL)

## example user-specified vaccine allocation function
# allocate vaccines according to absolute incidence in each location, then
# allocate uniformly within each country (without discriminating between
# ages/risk groups/infection statuses)
# needs to take the arguments
# sum_age_risk_func: a function which sums a vector across age and risk groups
# (created automatically in the code)
# travel_matrix: a square matrix with side length n_countries. 
# travel_matrix[x,y] is the proportion of time an individual in location x spends
# in location y
# vax_allocation_params: list of parameters for vaccine allocation
# S: numeric vector of length n. number of unvaccinated susceptibles
# E: numeric vector of length n. number of unvaccinated exposed
# I: numeric vector of length n. number of unvaccinated infectious
# R: numeric vector of length n. number of unvaccinated recovered
# vax_pool: numeric vector of length 1. number of vaccines available
user_specified_vax_alloc_func <- function(sum_age_risk_func, 
                                          travel_matrix,
                                          vax_allocation_params,
                                          S, E, I, R, vax_pool) {
  # find incidence in each country
  # incidence proportional to E
  incidence_by_country <- sum_age_risk_func(E)
  if(any(incidence_by_country > 0)) {
    n_vax_allocated <- incidence_by_country / sum(incidence_by_country) * vax_pool
  } else {
    n_vax_allocated <- incidence_by_country * 0 # allocate nothing
  }
  return(n_vax_allocated)
}

# example function of vaccine production:
# no vaccine produced until time vax_production_params[["detection_delay"]] + 
# vax_production_params[["production_delay"]], then
# constant production rate until max number of doses ever made reached,
# then no production
# needs to take the arguments
# vax_production_params: list of parameters for vaccine production
# t: scalar: time
user_specified_cum_vax_pool_func <- function(vax_production_params, t) {
  t_since_production <- t - (vax_production_params[["detection_delay"]] + 
                               vax_production_params[["production_delay"]])
  if(t_since_production < 0) {
    0
  } else {
    min(vax_production_params[["max_vax"]],
        t_since_production * vax_production_params[["production_rate"]])
  }
}

## SIMULATION OPTIONS
simulation_flags <- list(ageMixing=TRUE,
                         riskGroups=TRUE,
                         normaliseTravel=TRUE,
                         spatialCoupling=TRUE,
                         real_data = TRUE,
                         country_specific_contact = TRUE,
                         seasonal = FALSE,
                         rng_seed = 0)

## run simulation for tmax days
tmax <- 100
## tdiv timesteps per day
tdiv <- 24 ##DH

## Seasonality resolution
# 1 -> tmax*tdiv
# 12 -> tmax*tdiv/12
seasonality_resolution <- 1

## allocate and distribute vaccine every vac_alloc_period time divisions
## i.e. in this example, every 7 days
vax_alloc_period <- 24 * 7
n_riskgroups <- length(life_history_params$case_fatality_ratio)

# set random number generation seed
if(!is.null(simulation_flags[["rng_seed"]])) {
  set.seed(simulation_flags[["rng_seed"]])
}

if(simulation_flags[["real_data"]]) {
  ## get number of countries and ages from files
  demography_filename <- paste0(wd, "data/demographic_data_intersect.csv")
  contact_filename <- paste0(wd, "data/contact_data_intersect.csv")

  travel_filename <- paste0(wd, "data/flight_data_intersect.csv")
  latitude_filename <- paste0(wd, "data/latitudes_intersect.csv")
  risk_filename <- paste0(wd, "data/risk_group_data.csv")
    
  tmp <- read.csv(demography_filename, sep = ",")
  n_countries <- nrow(tmp)
  n_ages <- ncol(tmp) - 2
  risk_propns <- as.matrix(read.csv(risk_filename, sep = ",", header = FALSE))
  if(nrow(risk_propns) != n_ages) {
    stop("number of age groups inconsistent between data sets")
  }
  if(ncol(risk_propns) !=n_riskgroups) {
    stop("number of risk groups inconsistent between data sets")
  }
  
} else {
  ## SETUP FAKE COUNTRY DATA
  popn_size <- 100000
  ## Setup age propns
  n_ages <- 4
  age_propns <- rep(1/n_ages, n_ages)
  age_propns <- c(5,14,45,16)/80
  n_countries <- 10

  risk_propns <- rep(1/n_riskgroups,n_riskgroups) ## Assume risk groups are uniformly distributed
  risk_propns <- matrix(rep(risk_propns,each=n_ages),ncol=n_riskgroups) ## Assume that proportion of ages in each risk group are the same for all ages
}

## Setup risk groups

risk_factors <- rep(1, n_riskgroups) ## Assume that each risk group has same modifier

## Enumerate out risk factors for each age group
age_specific_riskgroup_factors <- matrix(rep(risk_factors,each=n_ages),
                                         ncol=n_riskgroups)

case_fatality_ratio_vec <- expand.grid("Location"=seq_len(n_countries), 
                                      "case_fatality_ratio" = life_history_params$case_fatality_ratio,
                                      "Age"=seq_len(n_ages))
case_fatality_ratio_vec <- case_fatality_ratio_vec$case_fatality_ratio

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

## process the vaccine production function
cum_vax_pool_func <- cum_vax_pool_func_closure(user_specified_cum_vax_pool_func, vax_production_params)

## process the vaccine allocation function

vax_allocation_func <- vaccine_allocation_closure(user_specified_vax_alloc_func, K, vax_allocation_params, labels)

## gather simulation parameters
sim_params <- list(n_countries=n_countries,
                   n_ages=n_ages,
                   n_riskgroups=n_riskgroups,
                   seed_vec = seed_vec,
                   seasonality_resolution=seasonality_resolution)

## run simulation
simulation_flags$seasonal <- TRUE
tmax <- 365
sim_params$seasonality_resolution <- tmax*tdiv/12
#sim_params$seasonality_resolution <- 1
Rprof(tmp <- tempfile(), line.profiling=TRUE)
system.time(
res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                      case_fatality_ratio_vec, X, labels, C3, K, latitudes, 
                      cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period)
)
Rprof()
summaryRprof(tmp, lines="show")
library(proftools)
plotProfileCallGraph(readProfileData(tmp),score = "total")
## plot stuff 
plot_labels <- expand.grid("Time"=seq(0,tmax,by=1/tdiv),"Location"=1:n_countries,"Age"=1:n_ages,"RiskGroup"=1:n_riskgroups)

tend <- ncol(res$S)
deaths <- X - res$S[,tend] - res$SV[,tend] - res$E[,tend] - res$EV[,tend] -
  res$I[,tend] - res$IV[,tend] - res$R[,tend] - res$RV[,tend]

popTotal <- sum(X)
globalAttack <- sum(res$R[,tend] + res$RV[,tend] + deaths)/popTotal

I <- cbind(labels[,c("Location","Age","RiskGroup")], res$I + res$IV)
I <- melt(I, id.vars=c("Location","Age","RiskGroup"))
I$Age <- as.factor(I$Age)
I$RiskGroup <- as.factor(I$RiskGroup)
I$variable <- as.numeric(I$variable)
times <- seq(0,tmax,by=1/tdiv)
I$variable <- times[I$variable]
#I_aggregated <- aggregate(I[,"value"], I[,c("variable","Location","Age")], FUN=sum)
I <- data.table(I)
setkey(I, "variable","Location","Age")
I_aggregated <- I[,x:=sum(value),by=key(I)]
N <- data.table(aggregate(data=labels, X~Location + Age,FUN=sum))
#N <- aggregate(data=labels, X~Location + Age,FUN=sum)
N <- N[N$Location %in% unique(I_aggregated$Location),]

I_aggregated <- merge(I_aggregated,N,by=c("Location","Age"))

p_normal <- ggplot(I_aggregated[I_aggregated$Location %in% unique(I_aggregated$Location)[1:20],],aes(x=variable,y=x/X,col=Age)) +
    geom_line() +
    facet_wrap(~Location) +
    theme_bw()
pdf("monthly.pdf")
plot(p_normal)
dev.off()
# p2 <- ggplot(I, aes(x=variable,y=value,col=RiskGroup)) + geom_line() + facet_grid(Age~Location) + theme_bw()

# to do: make plot of vaccines allocated (SV + EV + IV + RV) over time by country
# grid_plot <- cowplot::plot_grid(p1,p2,ncol=2,align="hv")
