#' Cluster submission for amp and epsilon calibration
#' 
#' Function used to submit jobs to the cluster, solving the global model with specified parameters. This particular function is intended to test different values of amp and epsilon, which is why they are amongst the first arguments. To write your own function which can be submitted to the cluster, simply place the parameters that you want to change as the first few arguments.
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
#' @param requested_stats character string to be passed to calculate_summaries function, to tell it which summary statistics to calculate and save
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
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
                                       seed_params, requested_stats, other_info){
  ## Set seed if specified
  if(!is.null(simulation_flags[["rng_seed"]])) {
    set.seed(simulation_flags[["rng_seed"]])
  }
  
  ## change seasonality amplitude and travel parameter
  travel_params[["epsilon"]] <- epsilon
  seasonality_params[["amp"]] <- amp
  
  #############################################################################
  # the below lines should be the same for all functions to run on the cluster
  #############################################################################
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
                        n_runs=n_runs, requested_stats=requested_stats, other_info)
  message("Simulations complete")
  
  #############################################################################
  # the above lines should be the same for all functions to run on the cluster
  #############################################################################
  # postprocessing: plot peak times
  labels <- processed_inputs[["labels"]]
  regionDat <- other_info[["regionDat"]]
  latitudeDat <- other_info[["latitudeDat"]]
  # plot_peak_times outputs a plot and a data frame containing the output
  final <- plot_peak_times(res, labels, 
                           regionDat, 
                           latitudeDat)
  p <- final[[1]]
  dat <- final[[2]]
  
  # save the plot
  filename <- paste0("outputs/calibration_",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  write.table(dat, paste0(filename, "_data.csv"),sep=",",row.names=FALSE)
  return(TRUE)    
}

#' Cluster submission with fixed parameters
#' 
#' Function used to submit jobs to the cluster, solving the global model with specified parameters. This particular function is intended to test different values of amp and epsilon, which is why they are amongst the first arguments. To write your own function which can be submitted to the cluster, simply place the parameters that you want to change as the first few arguments.
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
#' @param requested_stats character string to be passed to calculate_summaries function, to tell it which summary statistics to calculate and save
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @return returns TRUE if the routine runs correctly. Will create a .csv file of peak time summaries for this run, and a .png plotting the distribution of peak times by country. The output filenames are:
#' 1) outputs/default_data.csv
#' 2) outputs/default_plot.png
#' @export
run_fixed_params <- function(n_runs, time_params, seasonality_params,
                                       life_history_params, travel_params, simulation_flags,
                                       vax_params,vax_production_params, vax_allocation_params, 
                                       user_specified_cum_vax_pool_func,
                                       user_specified_vax_alloc_func,
                                       seed_params, requested_stats, other_info){
  ## Set seed if specified
  if(!is.null(simulation_flags[["rng_seed"]])) {
    set.seed(simulation_flags[["rng_seed"]])
  }
  
  #############################################################################
  # the below lines should be the same for all functions to run on the cluster
  #############################################################################
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
                        n_runs=n_runs, requested_stats=requested_stats, other_info)
  message("Simulations complete")
  
  #############################################################################
  # the above lines should be the same for all functions to run on the cluster
  #############################################################################
  # postprocessing: plot peak times
  labels <- processed_inputs[["labels"]]
  regionDat <- other_info[["regionDat"]]
  latitudeDat <- other_info[["latitudeDat"]]
  # plot_peak_times outputs a plot and a data frame containing the output
  final <- plot_peak_times(res, labels, 
                           regionDat, 
                           latitudeDat)
  p <- final[[1]]
  dat <- final[[2]]
  
  # save the plot
  filename <- "outputs/default"
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  write.table(dat, paste0(filename, "_data.csv"),sep=",",row.names=FALSE)
  return(TRUE)    
}