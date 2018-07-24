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

# load previous simulation arguments
previous_filename <- "raw_data/example_args_list.rds"
previous_args <- readRDS(previous_filename)
list2here(previous_args)

################################################################################
# change whatever parameters you want here

# Where to save simulation results
outputDir <- "outputs"
if(!file.exists(outputDir)) dir.create(outputDir)
output_prefix <- "repeats"
output_prefix <- paste(outputDir, output_prefix, sep = "/")

# character string specifying function used to calculate summaries of each run.
# see current options in get_vaxedemic_func_options()
# when writing these functions, the argument names must be things that can found in the environment
# after running the main simulation
calculate_summaries_func <- "calc_incidence_vaccinated_peak_times_attack_rates"

# character string specifying function to do postprocessing
# see current options in get_vaxedemic_func_options()
postprocessing_func <- "postprocessing_incidence_vaccinated_peak_times_attack_rates"

n_runs <- 3

################################################################################

# if TRUE, run a short test before the full number of runs
short_test <- FALSE
# if cluster && (!test_local), run test on cluster, otherwise run locally
test_local <- TRUE
# number of runs per test
n_runs_test <- 2
# if n_runs <= n_runs_test, don't run teh test regardless of teh value of short_test set above
short_test <- short_test && (n_runs > n_runs_test)

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
    saveRDS(args_list, paste0(output_prefix,"_args_list.rds"))
    
    if(short_test) {
      if(test_local) {
        args_list_temp <- make_arg_list(runs = NULL, run_func, obj = NULL)
        args_list_temp <- shorten_runs(args_list_temp, n_runs_test)
        do.call(run_func, args_list_temp)
      } else {
        args_list_temp <- args_list
        args_list_temp <- shorten_runs(args_list_temp, n_runs_test)
        job_test <- obj$enqueue(do.call(run_func, args_list_temp))
      }
    }

    job <- obj$enqueue(do.call(run_func, args_list))
  } else {
    # run a single job
    args_list <- make_arg_list(runs = NULL, run_func, obj = NULL)
    saveRDS(args_list, paste0(output_prefix,"_args_list.rds"))
    
    if(short_test) {
      args_list_temp <- args_list
      args_list_temp <- shorten_runs(args_list_temp, n_runs_test)
      do.call(run_func, args_list_temp)
    }
    
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
    saveRDS(args_list, paste0(output_prefix,"_args_list.rds"))
    
    if(short_test) {
      if(test_local) {
        # run test for first parameter set only
        args_list_temp <- make_arg_list(runs, run_func, obj = NULL)
        args_list_temp <- args_list_temp[[1]]
        args_list_temp <- shorten_runs(args_list_temp, n_runs_test)
        do.call(run_func, args_list_temp)
      } else {
        args_list_temp <- args_list[[1]]
        args_list_temp <- shorten_runs(args_list_temp, n_runs_test)
        job_test <- obj$enqueue(do.call(run_func, args_list_temp))
      }
    }
    
    jobs <- do.call(queuer::enqueue_bulk, args_list)
  } else {
    # run locally
    args_list <- make_arg_list(runs, run_func, obj = NULL)
    saveRDS(args_list, paste0(output_prefix,"_args_list.rds"))
    
    if(short_test) {
        args_list_temp <- args_list[[1]]
        args_list_temp <- shorten_runs(args_list_temp, n_runs_test)
        do.call(run_func, args_list_temp)
    }
    
    lapply(args_list, function(x) do.call(run_func, x))
  }
}
