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
calibrating_amp_and_travel <- function(amp, epsilon,regionDat, latitudeDat, simulation_flags, life_history_params, vax_params, sim_params,
                                       case_fatality_ratio_vec, popns, labels, contactMatrix, travelMatrix, latitudes, 
                                       user_specified_cum_vax_pool_func, vax_production_params,
                                       user_specified_vax_alloc_func, vax_allocation_params,
                                       tmax, tdiv, vax_alloc_period,
                                       n_runs=5){
    ## process the vaccine production function
    cum_vax_pool_func <- cum_vax_pool_func_closure(user_specified_cum_vax_pool_func, vax_production_params)

    ## process the vaccine allocation function
    vax_allocation_func <- vaccine_allocation_closure(user_specified_vax_alloc_func,
                                                      travelMatrix, vax_allocation_params, labels)
    
    
    sim_params[["amp"]] <- amp
    travel_params[["epsilon"]] <- epsilon    
    res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                          case_fatality_ratio_vec, popns, labels, contactMatrix, travelMatrix, latitudes, 
                          cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period,
                          n_runs=5)
    
    final <- plot_peak_times(res, labels, regionDat, latitudeDat)
    p <- final[[1]]
    dat <- final[[2]]

    filename <- paste0("outputs/calibration_",amp,"_",epsilon)
        
    png(paste0(filename, "_plot.png",width=800,height=1200))
    plot(p)
    dev.off()
    write.table(dat, paste0(filename, "_data.csv"),sep=",",row.names=FALSE)
    return(TRUE)    
}