library(reshape2)
library(ggplot2)
library(Matrix)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=4)
wd <- "~/Documents/vaxedemic/" 
devtools::load_all(wd)

###################################################
## RUN CONTROL
###################################################
## Number of simulations to run
n_runs <- 1
## run simulation for tmax days
tmax <- 365
## tdiv timesteps per day
tdiv <- 24
## Seasonality resolution
seasonality_resolution <- tmax*tdiv/12 # Average seasonality into 12 blocks of time
###################################################


###################################################
## LIFE HISTORY PARAMETER INPUTS
###################################################
## R_0, recovery time and latent period
life_history_params <- list(R0=1.4, TR=2.6, LP = 1.5, case_fatality_ratio = rep(2e-2,2))

## travel parameters: scaling of off-diagonals
travel_params <- list(epsilon = 1e-3)
###################################################


###################################################
## SIMULATION OPTIONS
###################################################
simulation_flags <- list(ageMixing=TRUE,
                         riskGroups=TRUE,
                         normaliseTravel=TRUE,
                         spatialCoupling=TRUE,
                         real_data = TRUE,
                         country_specific_contact = TRUE,
                         seasonal = TRUE,
                         rng_seed = 0)
###################################################


###################################################
## VACCINE INPUTS
###################################################
## vaccine efficacy and initial vaccinated proportion

vax_params <- list(efficacy = .7, propn_vax0 = 0)

## vaccine production function
user_specified_cum_vax_pool_func <- produce_vax_linear_with_delay
## parameters for vaccine production (see cum_vax_pool_func_closure)
vax_production_params <- list(detection_delay = 0, production_delay = 365/2, 
                              production_rate = 550e6/(365/12*3), max_vax = Inf)

## vaccination strategy
user_specified_vax_alloc_func <- vaccinate_by_incidence
## parameters for vaccine allocation
vax_allocation_params <- list(priorities = NULL)

## allocate and distribute vaccine every vac_alloc_period time divisions
## i.e. in this example, every 7 days
vax_alloc_period <- 24 * 7
###################################################


###################################################
## SETUP SEEDING
###################################################
# set random number generation seed
if(!is.null(simulation_flags[["rng_seed"]])) {
  set.seed(simulation_flags[["rng_seed"]])
}

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


###################################################
## SETUP INPUTS
###################################################
all_inputs <- setup_inputs(wd, simulation_flags, life_params ,travel_params)

## Extract inputs from setup
popns <- all_inputs$popns
labels <- all_inputs$labels
contactMatrix <- all_inputs$contactMatrix
travelMatrix <- all_inputs$travelMatrix
latitudes <- all_inputs$latitudes
n_countries <- all_inputs$n_countries
n_ages <- all_inputs$n_ages
n_riskgroups <- all_inputs$n_riskgroups
###################################################

###################################################
## FINAL SETUP
###################################################
case_fatality_ratio_vec <- expand.grid("Location"=seq_len(n_countries), 
                                      "case_fatality_ratio" = life_history_params$case_fatality_ratio,
                                      "Age"=seq_len(n_ages))
case_fatality_ratio_vec <- case_fatality_ratio_vec$case_fatality_ratio

#construct vector of number of exposed individuals in each location, age, risk group 
## for now, can only seed in one location, age, risk group. Vectorise later
seed_vec <- double(length(popns))
seed_vec[(which(labels$Location == seedCountries & labels$Age == seedAges &
              labels$RiskGroup == seedRiskGroups))[1]] <- seedSizes

## process the vaccine production function
cum_vax_pool_func <- cum_vax_pool_func_closure(user_specified_cum_vax_pool_func, vax_production_params)

## process the vaccine allocation function
vax_allocation_func <- vaccine_allocation_closure(user_specified_vax_alloc_func,
                                                  travelMatrix, vax_allocation_params, labels)

## gather simulation parameters
sim_params <- list(n_countries=n_countries,
                   n_ages=n_ages,
                   n_riskgroups=n_riskgroups,
                   seed_vec = seed_vec,
                   seasonality_resolution=seasonality_resolution,
                   tdelay=180,
                   amp=.7)
###################################################
res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                      case_fatality_ratio_vec, popns, labels, contactMatrix, travelMatrix, latitudes, 
                      cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period,
                      n_runs=5)
