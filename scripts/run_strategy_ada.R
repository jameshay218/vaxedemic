run_strategy <- function(strategy, production_delay, stockpile_size = 0, seedCountries = "China") {
  
  stopifnot(strategy %in% c("no_vaccination", 
                            "incidence", 
                            "curr_alloc",
                            "top_n_countries",
                            "fn_pop_size"))
  
  cluster <- TRUE # run on cluster or locally
  # user identifier -- only needed if running on cluster
  user <- "ayan"
  
  # if TRUE, run for one fixed set of parameters;
  # if FALSE, run for many combinations of parameters
  run_fixed <- strategy %in% c("no_vaccination", "incidence", "curr_alloc")
  
  # load vaxedemic package
  # local directory with the vaxedemic package
  package_dir <- "~/Documents/vaxedemic/"
  devtools::load_all(package_dir)
  setwd(package_dir)
  
  # Where to save simulation results
  outputDir <- paste0("temp_outputs2/pd",
                      production_delay, 
                      strategy, 
                      "_stockpile",
                      num2str(stockpile_size),
                      "_",
                      seedCountries)
  
  output_prefix <- strategy
  output_prefix <- paste(outputDir, output_prefix, sep = "/")
  
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
  n_runs <- 500
  
  # parameters to do with time steps in simulation
  time_params <- list(tmax = 720, # Maximum time of simulation (in days) -- 
                      # needs to be multiple of seasonality_params[["days_per_block"]]
                      tdiv = 6) # Number of time steps per day
  
  # seasonality parameters
  seasonality_params <- list(tdelay = 180, # (in days) Shifts the seasonality function - changing this effectively changes the seed time.
                             # tdelay = 0 is seed at t = 0 in sinusoidal curve, roughly start of autumn in Northern hemisphere
                             days_per_block = 30, # average seasonality over blocks of this many days
                             amp = 0.3) # amplitude of seasonality
  
  ## Life history parameters, including R0
  life_history_params <- list(R0=1.4, TR=2.6, LP = 1.5, case_fatality_ratio = rep(1e-3,2))
  
  ## Travel parameters
  travel_params <- list(epsilon = 5e-5)
  
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
  vax_params <- list(efficacy = .7)
  
  # parameters to do with vaccine production. correspond to arguments of user_specified_cum_vax_pool_func
  vax_production_params <- list(detection_delay = 0, production_delay = production_delay, 
                                production_rate = 550e06/(365/12*3), max_vax = Inf,
                                stockpile_size = stockpile_size)
  if(strategy == "no_vaccination") {
    vax_production_params$production_delay <- Inf
  }
  # parameters to do with vaccine allocation. correspond to arguments of user_specified_vax_alloc_func
  vax_allocation_params <- list(priorities = NULL, period = 6 * 7, coverage = NULL)
  
  # name of vaccine production function in vaxedemic package.  must specify as character string for do.call to work
  # see current options in get_vaxedemic_func_options()
  user_specified_cum_vax_pool_func <- "produce_vax_linear_with_delay"
  # name of vaccine allocation function in vaxedemic package.  must specify as character string for do.call to work
  # see current options in get_vaxedemic_func_options()
  if(strategy == "incidence") {
    user_specified_vax_alloc_func <- "vaccinate_by_incidence"
  } else {
    user_specified_vax_alloc_func <- "vaccinate_by_current_seasonal_alloc"
  }

  
  # parameters to do with seeding the pandemic
  seedCountries <- seedCountries

  seed_params <- list(Countries = seedCountries, # where to seed
                      Sizes = c(20), # how many to seed in each country
                      Ages = 3, # which age group to seed in each country
                      RiskGroups = 1) # which risk group to seed in each country
  
  # character string specifying function used to calculate summaries of each run.
  calculate_summaries_func <- "calc_incidence_peak_times_attack_rates_deaths"
  
  # character string specifying function to do postprocessing
  postprocessing_func <- "postprocessing_incidence_peak_times_attack_rates_deaths"
  
  # other_info provides any other information needed, such as to calculate the summaries
  # or post-process results.
  # in this case, we need region and latitude information for each country to make the plots.
  regionDat <- read.csv("data/regions_clean.csv")
  latitudeDat <- read.csv("data/latitudes_intersect.csv")
  other_info <- list(regionDat = regionDat,
                     latitudeDat = latitudeDat)
  
  if(cluster) {
    # Setup an interface to the cluster
    # sometimes fails with "Error in buildr_http_client_response(r) : Not Found (HTTP 404)" -- just re-run
    obj <- setup_cluster(user,expire=1e-10)
  } else if(.Platform$OS.type == "unix") {
    library(doMC)
    registerDoMC(cores=4)
  }
  
  if(!file.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

  if(run_fixed) {
    ################################################################################
    # run for fixed parameters
    ################################################################################
    run_func <- "run_fixed_params_and_postprocess"
    if(cluster) {
      # submit to cluster
      args_list <- make_arg_list(runs = NULL, run_func, obj)
      saveRDS(args_list, paste0(output_prefix,"_args_list.rds"))
      write(Sys.time(), paste0(output_prefix,"_timestamp.txt"))
      job <- obj$enqueue(do.call(run_func, args_list))
    } else {
      # run a single job
      args_list <- make_arg_list(runs = NULL, run_func, obj = NULL)
      saveRDS(args_list, paste0(output_prefix,"_args_list.rds"))
      write(Sys.time(), paste0(output_prefix,"_timestamp.txt"))
      
      do.call(run_func, args_list)
    }
  } else {
    ################################################################################
    # run for different combinations of parameters
    ################################################################################
    
    # the function to be run to vary parameters.
    run_func <- "random_vacc_allocations"
    
    ## Generate all combinations of these two parameters
    ## Generate a data frame for these parameters and a run name
    ## identifier for each combination. The column names for this 
    ## data frame must correspond to the first 
    ## arguments of run_func
    if(strategy == "top_n_countries") {
      data_dir <- "data/coverage_tables_by_popn/"
      n_countries <- c(1, 2, 5, 10, 127)
      filenames <- paste0("coverage_data_", n_countries, ".csv")
    } else {
      data_dir <- "data/coverage_tables_fn_pop_size/"
      alpha <- c(.5, 1, 2, 3)
      filenames <- paste0("coverage_data_", num2str(alpha), ".csv")
    }

    filenames_no_ext <- vcapply(filenames,
                                function(x) substr(x, 1, nchar(x) - 4))
    
    runs <- data.frame("runName"=filenames_no_ext,
                       "coverage_filename"=paste0(data_dir, filenames),
                       stringsAsFactors = FALSE)
    runs$runName <- as.character(runs$runName)
    
    # run in cluster or locally
    if(cluster) {
      # submit to cluster
      args_list <- make_arg_list(runs, run_func, obj)
      saveRDS(args_list, paste0(output_prefix,"_args_list.rds"))
      write(Sys.time(), paste0(output_prefix,"_timestamp.txt"))
      jobs <- do.call(queuer::enqueue_bulk, args_list)
    } else {
      # run locally
      
      args_list <- make_arg_list(runs, run_func, obj = NULL)
      saveRDS(args_list, paste0(output_prefix,"_args_list.rds"))
      write(Sys.time(), paste0(output_prefix,"_timestamp.txt"))
      lapply(args_list, function(x) do.call(run_func, x))
    }
  }
}

run_all_strategies <- function() {

  run_strategy("no_vaccination", 0, 0)
  strategy <- c("incidence", 
                "curr_alloc",
                "top_n_countries",
                "fn_pop_size")
  
  production_delay <- c(7, 90, 180)
  
  stockpile_size <- 550e6
  
  # seedCountries <- c("China", "Sao_Tome_and_Principe", "Belgium",
  #                    "Singapore", "Uganda")
  seedCountries <- c("Sao_Tome_and_Principe", "Belgium",
                     "Singapore", "Uganda")
  # largest pop size, smallest pop size, medium pop size, 
  # 3rd largest connnectivity, smallest connectivity
  
  pars <- expand.grid(strategy = strategy, production_delay = production_delay,
                      stringsAsFactors = FALSE)

  pars$stockpile_size <- 0
  pars$seedCountries <- "China"
  pars_stockpile <- data.frame(strategy = strategy, production_delay = 0, stockpile_size = stockpile_size)
  pars_stockpile <- expand.grid(strategy = c(strategy, "no_vaccination"), 
                                production_delay = 0, 
                                stockpile_size = stockpile_size,
                                seedCountries = seedCountries,
                                stringsAsFactors = FALSE)
  pars <- pars_stockpile
  # pars <- rbind(pars, pars_stockpile)
  Map(run_strategy, pars$strategy, pars$production_delay, pars$stockpile_size, pars$seedCountries)
}