library(reshape2)
library(ggplot2)
library(Matrix)
library(data.table)
wd <- "~/Documents/vaxedemic/" 
devtools::load_all(wd)

n_runs <- 5
tmax <- 365
tdiv <- 24
seasonality_resolution <- tmax*tdiv/12 # Average seasonality into 12 blocks of time
life_history_params <- list(R0=1.4, TR=2.6, LP = 1.5, case_fatality_ratio = rep(2e-2,2))
travel_params <- list(epsilon = 1e-3)
simulation_flags <- list(ageMixing=TRUE,
                         riskGroups=TRUE,
                         normaliseTravel=TRUE,
                         spatialCoupling=TRUE,
                         real_data = TRUE,
                         country_specific_contact = TRUE,
                         seasonal = TRUE,
                         rng_seed = 0)
vax_params <- list(efficacy = .7, propn_vax0 = 0)
vax_production_params <- list(detection_delay = 0, production_delay = 365/2, 
                              production_rate = 550e6/(365/12*3), max_vax = Inf)
vax_allocation_params <- list(priorities = NULL)
vax_alloc_period <- 24 * 7
if(simulation_flags[["real_data"]]) {
  seedCountries <- "China"
} else {
  seedCountries <- 1
}
seedSizes <- c(20)
seedAges <- 3
seedRiskGroups <- 1

regionDat <- read.csv("data/regions_clean.csv")
latitudeDat <- read.csv("data/latitudes_intersect.csv")

obj <- setup_cluster_JH("~/net/home/vaxedemic")

tdelay <- 180

epsilons <- c(0.00001, 0.00005, 0.0001,0.0005,0.001,0.005,0.01)
amps <- seq(0.1,1,by=0.05)
runs <- expand.grid(amp=amps, epsilon=epsilons)
runs <- cbind("runName"=paste0("test",1:nrow(runs)),runs)
runs$runName <- as.character(runs$runName)
jobs <- queuer::enqueue_bulk(obj, runs, "calibrating_amp_and_travel","", n_runs=100, tmax=tmax, tdiv=tdiv, 
                             seasonality_resolution=seasonality_resolution,
                             life_history_params=life_history_params, travel_params=travel_params, simulation_flags=simulation_flags,
                             vax_params=vax_params,vax_production_params=vax_production_params, vax_allocation_params=vax_allocation_params, 
                             vax_alloc_period=vax_alloc_period,
                             seedCountries=seedCountries, seedSizes=seedSizes, seedAges=seedAges, seedRiskGroups=seedRiskGroups,
                             tdelay=tdelay, regionDat=regionDat, latitudeDat=latitudeDat,do_call=TRUE,timeout=0)