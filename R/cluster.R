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





#' @export
calibrating_amp_and_travel <- function(runName, amp, epsilon, wd,
                                       n_runs, tmax, tdiv, seasonality_resolution,
                                       life_history_params, travel_params, simulation_flags,
                                       vax_params,vax_production_params, vax_allocation_params, vax_alloc_period,
                                       seedCountries, seedSizes, seedAges, seedRiskGroups,
                                       tdelay, regionDat, latitudeDat, requested_stats){
    if(!is.null(simulation_flags[["rng_seed"]])) {
        set.seed(simulation_flags[["rng_seed"]])
    }

    travel_params[["epsilon"]] <- epsilon
    
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
    
    final <- plot_peak_times(res[[1]], labels, regionDat, latitudeDat)
    p <- final[[1]]
    dat <- final[[2]]

    filename <- paste0("outputs/calibration_",runName)

    png(paste0(filename, "_plot.png"),width=800,height=1200)
    plot(p)
    dev.off()
    
    write.table(dat, paste0(filename, "_data.csv"),sep=",",row.names=FALSE)
    return(TRUE)    
}
