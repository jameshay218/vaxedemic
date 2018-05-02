#' Cluster setup for JH
#'
#' This function sets up an object linking to the DIDE cluster in a very crude way. Each user will need to implement their own version of this that returns a valid didehpc object to submit jobs to.
#' @param wd the working directory to run the cluster job from. This should be the user's network home drive eg. "~/net/home/vaxedemic"
#' @return a didehpc::queue_didehpc object
#' @export
setup_cluster_JH <- function(wd){
    setwd(wd)
    options(didehpc.credentials = "~/.smbcredentials",
            didehpc.cluster = "fi--didemrchnb")
    src <- provisionr::package_sources(local = "~/Documents/vaxedemic/")
    sources <- NULL
    
    ## Setup contexts
    context::context_log_start()
    root <- "contexts"
    packages <- list(attached=c("vaxedemic","plyr","reshape2","data.table","ggplot2","Matrix","foreach"))
    ctx <- context::context_save(packages=packages,path=root, sources=sources,package_sources=src)
    
    ## Submit setup to cluster
    obj1 <- didehpc::queue_didehpc(ctx)
    return(obj1)
}



#' Cluster submission for amp and epsilon calibration
#'
#' Function used to submit jobs to the cluster, solving the global model with specified parameters. This particular function is intended to test different values of amp and epsilon, which is why they are amongst the first arguments. To write your own function which can be submitted to the cluster, simply place the parameters that you want to change as the first few arguments.
#' @param runName the base of the filename to save outputs to. Will be appended to "outputs/calibration_". DO NOT SPECIFY A FILE EXTENSION - THIS IS HANDLED IN THE FUNCTION
#' @param amp the amplitude of seasonal forcing
#' @param epsilon connectivity of off-diagonal elements in the travel matrix
#' @param wd the working directory to run the cluster job from. This should be available to the DIDE network ie. your network Q: drive. Note that this should have all of the files and subfolders needed to run \code{\link{setup_inputs}}
#' @param n_runs number of simulations to run for this set of parameters
#' @param tmax the number of days for which to run the simulation
#' @param tdiv the number of timesteps per day
#' @param seasonality_resolution averages the seasonality function into buckets of this many time steps. For example, if you'd like 12 blocks of seasonality in this time frame, you should set this parameter to (tmax*tdiv)/12
#' @param life_history_params named vector (or list) with the numeric elements
#' "R0", "TR" (time to recovery) and "LP" (latent period)
#' @param travel_params list of parameters relating to travel. This should have the following elements: 1) epsilon, which scales the off-diagonals of the travel matrix
#' @param simulation_flags named vector (or list) with the logical elements "normaliseTravel"
#' and "seasonal"
#' @param vax_params named vector (or list) with the numeric elements "efficacy"
#' and "propn_vax0" (initial proportion of vaccinated individuals; assumed constant
#' across location, age and risk groups)
#' @param vax_production_params
#' @param vax_allocation_params
#' @param vax_allocation_period
#' @param seedCountries
#' @param seedSizes
#' @param seedAges
#' @param seedRiskGroups
#' @param tdelay
#' @param regionDat
#' @param latitudeDat
#' @param requested_stats
#' @return returns TRUE if the routine runs correctly. Will create a .csv file of peak time summaries for this run, and a .png plotting the distribution of peak times by country. The output filenames are:
#' 1) outputs/calibration_***runName***_data.csv
#' 2) outputs/calibration_***runName***_plot.png
#' @export
calibrating_amp_and_travel <- function(runName, amp, epsilon, wd,
                                       n_runs, tmax, tdiv, seasonality_resolution,
                                       life_history_params, travel_params, simulation_flags,
                                       vax_params,vax_production_params, vax_allocation_params, vax_alloc_period,
                                       seedCountries, seedSizes, seedAges, seedRiskGroups,
                                       tdelay, regionDat, latitudeDat, requested_stats){
    ## Set seed if specified
    if(!is.null(simulation_flags[["rng_seed"]])) {
        set.seed(simulation_flags[["rng_seed"]])
    }

    ## Add the given epsilon to the travel parameter list
    travel_params[["epsilon"]] <- epsilon

    ## Generate 
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

    seed_vec <- double(length(popns))
    seed_vec[(which(labels$Location == seedCountries & labels$Age == seedAges &
                    labels$RiskGroup == seedRiskGroups))[1]] <- seedSizes

    ## process the vaccine production function
    cum_vax_pool_func <- cum_vax_pool_func_closure(user_specified_cum_vax_pool_func, vax_production_params)

    ## process the vaccine allocation function
    vax_allocation_func <- vaccine_allocation_closure(user_specified_vax_alloc_func,
                                                      travelMatrix, vax_allocation_params, labels)

    sim_params <- list(n_countries=n_countries,
                       n_ages=n_ages,
                       n_riskgroups=n_riskgroups,
                       seed_vec = seed_vec,
                       seasonality_resolution=seasonality_resolution,
                       tdelay=tdelay,
                       amp=amp)
    
    message("Setup complete")
    message(cat("Number of runs: ", n_runs,sep="\t"))
    
    res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                          case_fatality_ratio_vec, popns, labels, contactMatrix, travelMatrix, latitudes, 
                          cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period,
                          n_runs=n_runs, requested_stats=requested_stats, regionDat=regionDat, latitudeDat=latitudeDat)
    message("Simulations complete")

    ## Cbind country incidence
    ## 
    
    final <- plot_peak_times(res, labels, regionDat, latitudeDat)
    p <- final[[1]]
    dat <- final[[2]]

    filename <- paste0("outputs/calibration_",runName)

    png(paste0(filename, "_plot.png"),width=800,height=1200)
    plot(p)
    dev.off()
    
    write.table(dat, paste0(filename, "_data.csv"),sep=",",row.names=FALSE)
    return(TRUE)    
}
