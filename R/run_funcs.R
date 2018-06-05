#' submit job for amp and epsilon calibration
#' 
#' Function used to submit jobs, solving the global model with specified parameters. This particular function is intended to test different values of amp and epsilon, which is why they are amongst the first arguments. To write your own function which can be submitted to the cluster, simply place the parameters that you want to change as the first few arguments.
#' @param runName the base of the filename to save outputs to. Will be appended to "outputs/calibration_". DO NOT SPECIFY A FILE EXTENSION - THIS IS HANDLED IN THE FUNCTION
#' @param amp the amplitude of seasonal forcing
#' @param epsilon connectivity of off-diagonal elements in the travel matrix
#' @param n_runs number of simulations to run for this set of parameters
#' @param time_params list of parameters to do with time steps in simulation.
# contains the elements tmax (Maximum time of simulation), tdiv (Number of time steps per day)
#' @param seasonality_params list of seasonality parameters.
#' contains the elements tdelay (0 <= tdelay <= 364): shifts the seasonality function - changing this effectively changes the seed time.
#' tdelay = 0 is seed at t = 0 in sinusoidal curve, roughly start of autumn in Northern hemisphere
#' division: Average seasonality into this many blocks of time
#' amp: amplitude of seasonality
#' @param life_history_params named vector (or list) with the numeric elements
#' "R0", "TR" (time to recovery) and "LP" (latent period)
#' @param travel_params list of parameters relating to travel. This should have the following elements: 1) epsilon, which scales the off-diagonals of the travel matrix
#' @param simulation_flags named vector (or list) with the logical elements "normaliseTravel"
#' and "seasonal"
#' @param vax_params named vector (or list) with the numeric elements "efficacy"
#' and "propn_vax0" (initial proportion of vaccinated individuals; assumed constant
#' across location, age and risk groups)
#' @param vax_production_params named vector (or list) with named elements matching the arguments of user_specified_cum_vax_pool_func
#' @param vax_allocation_params named vector (or list) with named elements matching the arguments of user_specified_vax_alloc_func
#' @param user_specified_cum_vax_pool_func function or character string of function to produce vaccine pool
#' @param user_specified_vax_alloc_func function or character string of function to allocate vaccines
#' @param seed_params named vector (or list) of parameters to do with seeding the pandemic.
#' contains the elements seedCountries: vector of country names in which to seed
# Sizes: vector of how many to seed in each country
# Ages: vector of which age group to seed in each country
# RiskGroups: vector of which risk group to seed in each country
#' @param calculate_summaries_func string specifying what function to use to calculate summaries.
#' @param postprocessing_func character string specifying function to do postprocessing
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @param output_prefix character vector of length 1.  Prefix for output filenames
#' @return returns TRUE if the routine runs correctly. Will create a .csv file of peak time summaries for this run, and a .png plotting the distribution of peak times by country. The output filenames are:
#' 1) outputs/calibration_***runName***_data.csv
#' 2) outputs/calibration_***runName***_plot.png
#' @export
calibrating_amp_and_travel <- function(runName, amp, epsilon, 
                                       n_runs, time_params, seasonality_params,
                                       life_history_params, travel_params, simulation_flags,
                                       vax_params,vax_production_params, vax_allocation_params, 
                                       user_specified_cum_vax_pool_func,
                                       user_specified_vax_alloc_func,
                                       seed_params, calculate_summaries_func, 
                                       postprocessing_func, other_info,
                                       output_prefix){
  
  ## change seasonality amplitude and travel parameter
  travel_params[["epsilon"]] <- epsilon
  seasonality_params[["amp"]] <- amp
  
  # run simulations
  # the following two lines get the arguments of run_fixed_params from the environment
  # and runs run_fixed_params with those arguments
  args_list <- list_vars_from_environment(formalArgs("run_fixed_params"))
  res_list <- do.call("run_fixed_params", args_list)
  
  # postprocessing
  # the following two lines get the arguments of postprocessing_func from the environment
  # and runs postprocessing_func with those arguments
  args_list <- list_vars_from_environment(formalArgs(postprocessing_func))
  do.call(postprocessing_func, args_list)

  return(TRUE)    
}

#' submit job with fixed parameters and postprocessing
#' 
#' Function used to submit jobs, solving the global model with specified parameters. This particular function is intended to test different values of amp and epsilon, which is why they are amongst the first arguments. To write your own function which can be submitted to the cluster, simply place the parameters that you want to change as the first few arguments.
#' @param n_runs number of simulations to run for this set of parameters
#' @param time_params list of parameters to do with time steps in simulation.
# contains the elements tmax (Maximum time of simulation), tdiv (Number of time steps per day)
#' @param seasonality_params list of seasonality parameters.
#' contains the elements tdelay (0 <= tdelay <= 364): shifts the seasonality function - changing this effectively changes the seed time.
#' tdelay = 0 is seed at t = 0 in sinusoidal curve, roughly start of autumn in Northern hemisphere
#' division: Average seasonality into this many blocks of time
#' amp: amplitude of seasonality
#' @param life_history_params named vector (or list) with the numeric elements
#' "R0", "TR" (time to recovery) and "LP" (latent period)
#' @param travel_params list of parameters relating to travel. This should have the following elements: 1) epsilon, which scales the off-diagonals of the travel matrix
#' @param simulation_flags named vector (or list) with the logical elements "normaliseTravel"
#' and "seasonal"
#' @param vax_params named vector (or list) with the numeric elements "efficacy"
#' and "propn_vax0" (initial proportion of vaccinated individuals; assumed constant
#' across location, age and risk groups)
#' @param vax_production_params named vector (or list) with named elements matching the arguments of user_specified_cum_vax_pool_func
#' @param vax_allocation_params named vector (or list) with named elements matching the arguments of user_specified_vax_alloc_func
#' @param user_specified_cum_vax_pool_func function or character string of function to produce vaccine pool
#' @param user_specified_vax_alloc_func function or character string of function to allocate vaccines
#' @param seed_params named vector (or list) of parameters to do with seeding the pandemic.
#' contains the elements seedCountries: vector of country names in which to seed
# Sizes: vector of how many to seed in each country
# Ages: vector of which age group to seed in each country
# RiskGroups: vector of which risk group to seed in each country
#' @param calculate_summaries_func string specifying what function to use to calculate summaries.
#' @param postprocessing_func character string specifying function to do postprocessing
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @param output_prefix character vector of length 1.  Prefix for output filenames
#' @return returns TRUE if the routine runs correctly. Will create a .csv file of peak time summaries for this run, and a .png plotting the distribution of peak times by country. The output filenames are:
#' 1) outputs/default_data.csv
#' 2) outputs/default_plot.png
#' @export
run_fixed_params_and_postprocess <- function(n_runs, time_params, seasonality_params,
                                       life_history_params, travel_params, simulation_flags,
                                       vax_params,vax_production_params, vax_allocation_params, 
                                       user_specified_cum_vax_pool_func,
                                       user_specified_vax_alloc_func,
                                       seed_params, calculate_summaries_func, postprocessing_func,
                                       other_info, output_prefix) {
  
  runName <- "fixed" # need to provide a runName to make filenames for postprocessing
  
  # run simulations
  # the following two lines get the arguments of run_fixed_params from the environment
  # and runs run_fixed_params with those arguments
  args_list <- list_vars_from_environment(formalArgs("run_fixed_params"))
  res_list <- do.call("run_fixed_params", args_list)
  
  # postprocessing
  # the following two lines get the arguments of postprocessing_func from the environment
  # and runs postprocessing_func with those arguments
  args_list <- list_vars_from_environment(formalArgs(postprocessing_func))
  do.call(postprocessing_func, args_list)
  
  return(TRUE)    
}

#' submit job with fixed parameters
#' 
#' Function used to submit jobs, solving the global model with specified parameters. This particular function is intended to test different values of amp and epsilon, which is why they are amongst the first arguments. To write your own function which can be submitted to the cluster, simply place the parameters that you want to change as the first few arguments.
#' @param n_runs number of simulations to run for this set of parameters
#' @param time_params list of parameters to do with time steps in simulation.
# contains the elements tmax (Maximum time of simulation), tdiv (Number of time steps per day)
#' @param seasonality_params list of seasonality parameters.
#' contains the elements tdelay (0 <= tdelay <= 364): shifts the seasonality function - changing this effectively changes the seed time.
#' tdelay = 0 is seed at t = 0 in sinusoidal curve, roughly start of autumn in Northern hemisphere
#' division: Average seasonality into this many blocks of time
#' amp: amplitude of seasonality
#' @param life_history_params named vector (or list) with the numeric elements
#' "R0", "TR" (time to recovery) and "LP" (latent period)
#' @param travel_params list of parameters relating to travel. This should have the following elements: 1) epsilon, which scales the off-diagonals of the travel matrix
#' @param simulation_flags named vector (or list) with the logical elements "normaliseTravel"
#' and "seasonal"
#' @param vax_params named vector (or list) with the numeric elements "efficacy"
#' and "propn_vax0" (initial proportion of vaccinated individuals; assumed constant
#' across location, age and risk groups)
#' @param vax_production_params named vector (or list) with named elements matching the arguments of user_specified_cum_vax_pool_func
#' @param vax_allocation_params named vector (or list) with named elements matching the arguments of user_specified_vax_alloc_func
#' @param user_specified_cum_vax_pool_func function or character string of function to produce vaccine pool
#' @param user_specified_vax_alloc_func function or character string of function to allocate vaccines
#' @param seed_params named vector (or list) of parameters to do with seeding the pandemic.
#' contains the elements seedCountries: vector of country names in which to seed
# Sizes: vector of how many to seed in each country
# Ages: vector of which age group to seed in each country
# RiskGroups: vector of which risk group to seed in each country
#' @param calculate_summaries_func string specifying what function to use to calculate summaries.
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @return returns a list containing the elements
#' res: summary statistics for each run
#' processed_inputs: results of pre-processing which might be needed for downstream analysis
#' @export
run_fixed_params <- function(n_runs, time_params, seasonality_params,
                             life_history_params, travel_params, simulation_flags,
                             vax_params,vax_production_params, vax_allocation_params, 
                             user_specified_cum_vax_pool_func,
                             user_specified_vax_alloc_func,
                             seed_params, calculate_summaries_func, other_info) {
  ## Set seed if specified
  if(!is.null(simulation_flags[["rng_seed"]])) {
    set.seed(simulation_flags[["rng_seed"]])
  }
  ## setup inputs 
  processed_inputs <- setup_inputs(simulation_flags, 
                                   life_history_params, 
                                   travel_params, 
                                   seed_params, 
                                   user_specified_cum_vax_pool_func, 
                                   vax_production_params,
                                   user_specified_vax_alloc_func, 
                                   vax_allocation_params)
  
  message("Setup complete")
  message(cat("Number of runs: ", n_runs,sep="\t"))
  
  res <- run_simulation(simulation_flags, life_history_params, vax_params, seasonality_params,
                        time_params, vax_allocation_params[["period"]], processed_inputs,
                        n_runs, calculate_summaries_func, other_info)
  message("Simulations complete")
  return(list(res = res, processed_inputs = processed_inputs))
}