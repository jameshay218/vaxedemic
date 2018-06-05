cluster <- FALSE # run on cluster or locally
# user identifier -- only needed if running on cluster
user <- "ayan"

# if TRUE, run for one fixed set of parameters;
# if FALSE, run for many combinations of parameters
run_fixed <- TRUE

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
# seed_params, calculate_summaries_func, postprocessing_func, other_info

## How many runs for each set of simulations?
n_runs <- 2

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
                         rng_seed = NULL)

# parameters to do with properties of the vaccine: efficacy and initial number vaccinated
vax_params <- list(efficacy = .7, propn_vax0 = 0)
# parameters to do with vaccine production. correspond to arguments of user_specified_cum_vax_pool_func
vax_production_params <- list(detection_delay = 0, production_delay = 365/2, 
                              production_rate = 550e6/(365/12*3), max_vax = Inf)
# parameters to do with vaccine allocation. correspond to arguments of user_specified_vax_alloc_func
vax_allocation_params <- list(priorities = NULL, period = 24 * 7, coverage = NULL)

# name of vaccine production function in vaxedemic package.  must specify as character string for do.call to work
# see current options in get_vaxedemic_func_options()
user_specified_cum_vax_pool_func <- "produce_vax_linear_with_delay"
# name of vaccine allocation function in vaxedemic package.  must specify as character string for do.call to work
# see current options in get_vaxedemic_func_options()
user_specified_vax_alloc_func <- "vaccinate_by_current_seasonal_alloc"

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

# character string specifying function used to calculate summaries of each run.
# see current options in get_vaxedemic_func_options()
# when writing these functions, the argument names must be things that can found in the environment
# after running the main simulation
calculate_summaries_func <- "calc_peak_times_and_attack_rates"

# character string specifying function to do postprocessing
# see current options in get_vaxedemic_func_options()
postprocessing_func <- "postprocessing_peak_times_and_attack_rates"

# certain postprocessing funcs require certain summaries to be calculated -- check
if((postprocessing_func == "postprocessing_country_attack" && calculate_summaries_func != "calc_country_attack") ||
   (postprocessing_func == "postprocessing_plot_peak_times" && calculate_summaries_func != "calc_peak_times") ||
   (postprocessing_func == "postprocessing_peak_times_and_attack_rates" && calculate_summaries_func != "calc_peak_times_and_attack_rates")) {
  stop("postprocessing function does not match summary function")
}

# other_info provides any other information needed, such as to calculate the summaries
# or post-process results.
# in this case, we need region and latitude information for each country to make the plots.
regionDat <- read.csv("data/regions_clean.csv")
latitudeDat <- read.csv("data/latitudes_intersect.csv")
other_info <- list(regionDat = regionDat,
                   latitudeDat = latitudeDat)

output_prefix <- "my_prefix"

if(cluster) {
  # Setup an interface to the cluster
  # sometimes fails with "Error in buildr_http_client_response(r) : Not Found (HTTP 404)" -- just re-run
  obj <- setup_cluster(user)
} else if(.Platform$OS.type == "unix") {
  library(doMC)
  registerDoMC(cores=4)
}

if(run_fixed) {
  ################################################################################
  # run for fixed parameters
  ################################################################################
  run_func <- "run_fixed_params_and_postprocess"
  if(cluster) {
    # submit to cluster
    args_list <- make_arg_list(runs = NULL, run_func, obj)
    job <- obj$enqueue(do.call(run_func, args_list))
  } else {
    # run a single job
    args_list <- make_arg_list(runs = NULL, run_func, obj = NULL)
    do.call(run_func, args_list)
  }
} else {
  ################################################################################
  # run for different combinations of parameters
  ################################################################################
  
  # the function to be run to vary parameters. write your own in run_funcs.R
  # must specify as character string for do.call to work
  # see current options in get_vaxedemic_func_options()
  run_func <- "calibrating_amp_and_travel"
  
  # set up the variable parameters.
  # in this case, we change the travel connectivity and seasonality amplitude.
  epsilons <- c(0.00001, 0.00005)
  amps <- c(0.1,1)
  
  ## Generate all combinations of these two parameters
  ## Generate a data frame for these parameters and a run name
  ## identifier for each combination. The column names for this 
  ## data frame must correspond to the first 
  ## arguments of run_func
  runs <- expand.grid(amp=amps, epsilon=epsilons)
  runs <- cbind("runName"=paste0("test",1:nrow(runs)),runs)
  runs$runName <- as.character(runs$runName)
  
  # run in cluster or locally
  if(cluster) {
    # submit to cluster
    args_list <- make_arg_list(runs, run_func, obj)
    jobs <- do.call(queuer::enqueue_bulk, args_list)
  } else {
    # run a single job
    args_list <- make_arg_list(runs, run_func, obj = NULL)
    lapply(args_list, function(x) do.call(run_func, x))
  }
}