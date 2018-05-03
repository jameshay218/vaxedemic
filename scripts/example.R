library(reshape2)
library(ggplot2)
library(Matrix)
library(data.table)

cluster <- TRUE # run on cluster or locally
run_fixed <- FALSE # run for one fixed set of parameters, or many combinations of parameters

# user identifier -- only needed if running on cluster
user <- "ayan"
# load vaxedemic package
# local directory with the vaxedemic package
package_dir <- "~/Documents/vaxedemic/"
devtools::load_all(package_dir)
setwd(package_dir)

# set up the arguments to be passed to the function to be run, which
# are constant across the sets of simulations.
# these arguments are
# n_runs, time_params, seasonality_params,
# life_history_params, travel_params, simulation_flags,
# vax_params,vax_production_params, vax_allocation_params,
# user_specified_cum_vax_pool_func,
# user_specified_vax_alloc_func,
# seed_params, requested_stats, other_info

## How many runs for each set of simulations?
n_runs <- 3

# parameters to do with time steps in simulation
time_params <- list(tmax = 1000, # Maximum time of simulation
                    tdiv = 24) # Number of time steps per day

# seasonality parameters
seasonality_params <- list(tdelay = 180, # Shifts the seasonality function - changing this effectively changes the seed time.
                           # tdelay = 0 is seed at t = 0 in sinusoidal curve, roughly start of autumn in Northern hemisphere
                           # Average seasonality into 12 blocks of time
                           division = 12,
                           amp = 0.1) # amplitude of seasonality

## Life history parameters, including R0
life_history_params <- list(R0=1.4, TR=2.6, LP = 1.5, case_fatality_ratio = rep(2e-2,2))

## Travel parameters
travel_params <- list(epsilon = 1e-3)

## Options for the simulation
simulation_flags <- list(ageMixing=TRUE,
                         riskGroups=TRUE,
                         normaliseTravel=TRUE,
                         spatialCoupling=TRUE,
                         real_data = TRUE,
                         country_specific_contact = TRUE,
                         seasonal = TRUE,
                         rng_seed = 0)

# parameters to do with properties of the vaccine: efficacy and initial number vaccinated
vax_params <- list(efficacy = .7, propn_vax0 = 0)
# parameters to do with vaccine production. correspond to arguments of user_specified_cum_vax_pool_func
vax_production_params <- list(detection_delay = 0, production_delay = 365/2, 
                              production_rate = 550e6/(365/12*3), max_vax = Inf)
# parameters to do with vaccine allocation. correspond to arguments of user_specified_vax_alloc_func
vax_allocation_params <- list(priorities = NULL, period = 24 * 7)

# name of vaccine production function in vaxedemic package.  must specify as character string for do.call to work
user_specified_cum_vax_pool_func <- "produce_vax_linear_with_delay"
# name of vaccine allocation function in vaxedemic package.  must specify as character string for do.call to work
user_specified_vax_alloc_func <- "vaccinate_by_incidence"

# parameters to do with seeding the pandemic
if(simulation_flags[["real_data"]]) {
  seedCountries <- "China"
} else {
  seedCountries <- 1
}

seed_params <- list(Countries = seedCountries, # where to seed
                    Sizes = c(20), # how many to seed in each country
                    Ages = 3, # which age group to seed in each country
                    RiskGroups = 1) # which risk group to seed in each country

# parameter to be passed to calculate_summaries function, to tell it which summary statistics to calculate and save
requested_stats <- "peak_times"
# other_info provides any other information needed, such as to calculate the summaries
# or post-process results.
# in this case, we need region and latitude information for each country to make the plots.
regionDat <- read.csv("data/regions_clean.csv")
latitudeDat <- read.csv("data/latitudes_intersect.csv")
other_info <- list(regionDat = regionDat,
                   latitudeDat = latitudeDat)

if(run_fixed) {
  ################################################################################
  # run for fixed parameters
  ################################################################################
  submit_fn <- "run_fixed_params"
  if(cluster) {
    # Setup an interface to the cluster
    # sometimes fails with "Error in buildr_http_client_response(r) : Not Found (HTTP 404)" -- just re-run
    obj <- setup_cluster(user)
    # submit to cluster
    args_list <- make_arg_list(runs = NULL, submit_fn, obj)
    job <- obj$enqueue(do.call(submit_fn, args_list))
  } else {
    # run a single job
    library(doMC)
    registerDoMC(cores=4)
    args_list <- make_arg_list(runs = NULL, submit_fn, obj = NULL)
    do.call(submit_fn, args_list)
  }
} else {
  ################################################################################
  # run for different combinations of parameters
  ################################################################################
  
  # the function to be run to vary parameters. write your own in funcs_to_run_on_cluster.R.
  # must specify as character string for do.call to work
  submit_fn <- "calibrating_amp_and_travel"
  
  # set up the variable parameters.
  # in this case, we change the travel connectivity and seasonality amplitude.
  epsilons <- c(0.00001, 0.00005)#, 0.0001,0.0005,0.001,0.005,0.01)
  amps <- seq(0.1,1)#,by=0.05)
  
  ## Generate all combinations of these two parameters
  ## Generate a data frame for these parameters and a run name
  ## identifier for each combination. The column names for this 
  ## data frame must correspond to the first 
  ## arguments of submit_fn
  runs <- expand.grid(amp=amps, epsilon=epsilons)
  runs <- cbind("runName"=paste0("test",1:nrow(runs)),runs)
  runs$runName <- as.character(runs$runName)
  
  # run in cluster or locally
  if(cluster) {
    # Setup an interface to the cluster
    # sometimes fails with "Error in buildr_http_client_response(r) : Not Found (HTTP 404)" -- just re-run
    obj <- setup_cluster(user)
    # submit to cluster
    args_list <- make_arg_list(runs, submit_fn, obj)
    jobs <- do.call(queuer::enqueue_bulk, args_list)
  } else {
    # run a single job
    library(doMC)
    registerDoMC(cores=4)
    args_list <- make_arg_list(runs, submit_fn, obj = NULL)
    lapply(args_list, function(x) do.call(submit_fn, x))
  }
}




