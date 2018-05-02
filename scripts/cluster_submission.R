library(reshape2)
library(ggplot2)
library(Matrix)
library(data.table)

## Working directory with the vaxedemic package
wd <- "~/Documents/vaxedemic/" 
devtools::load_all(wd)
setwd(wd)

## How many runs for each set of simulations?
n_runs <- 100
tmax <- 1000 # Maximum time of simulation
tdiv <- 24 # Number of time steps per day
seasonality_resolution <- tmax*tdiv/12 # Average seasonality into 12 blocks of time

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

## Vaccine function parameters
vax_params <- list(efficacy = .7, propn_vax0 = 0)
vax_production_params <- list(detection_delay = 0, production_delay = 365/2, 
                              production_rate = 550e6/(365/12*3), max_vax = Inf)
vax_allocation_params <- list(priorities = NULL)
vax_alloc_period <- 24 * 7

## Where to seed the pandemic
if(simulation_flags[["real_data"]]) {
  seedCountries <- "China"
} else {
  seedCountries <- 1
}
seedSizes <- c(20)
seedAges <- 3
seedRiskGroups <- 1

## This files should be provided in the vaxedemic package. They are used to provide
## region and latitude information for each country
regionDat <- read.csv("data/regions_clean.csv")
latitudeDat <- read.csv("data/latitudes_intersect.csv")


## Setup an interface to the cluster
obj <- setup_cluster_JH("~/net/home/vaxedemic")

## Shifts the seasonality function - changing this effectively changes the seed time
tdelay <- 180

####################
## SEASONALITY CALIBRATION
## Values of travel connectivity and seasonality amplitude to test
epsilons <- c(0.00001, 0.00005, 0.0001,0.0005,0.001,0.005,0.01)
amps <- seq(0.1,1,by=0.05)

## Generate all combinations of these two parameters
## Generate a data frame for these parameters and a run name
## identifier for each combination. The column names for this 
## data frame must correspond to the first three
## arguments of "calibrating_amp_and_travel"
runs <- expand.grid(amp=amps, epsilon=epsilons)
runs <- cbind("runName"=paste0("test",1:nrow(runs)),runs)
runs$runName <- as.character(runs$runName)

## Submit jobs to the cluster. The function being run is "calibrating_amp_and_travel",
## where "runs" specifies the inputs for the first three arguments. The rest of the
## arguments should be named and correspond to the rest of the arguments in 
## calibrating_amp_and_travel
jobs <- queuer::enqueue_bulk(obj, runs, "calibrating_amp_and_travel","", n_runs=n_runs, tmax=tmax, tdiv=tdiv, 
                             seasonality_resolution=seasonality_resolution,
                             life_history_params=life_history_params, travel_params=travel_params, simulation_flags=simulation_flags,
                             vax_params=vax_params,vax_production_params=vax_production_params, vax_allocation_params=vax_allocation_params, 
                             vax_alloc_period=vax_alloc_period,
                             seedCountries=seedCountries, seedSizes=seedSizes, seedAges=seedAges, seedRiskGroups=seedRiskGroups,
                             tdelay=tdelay, regionDat=regionDat, latitudeDat=latitudeDat, requested_stats="peak_times", do_call=TRUE,timeout=0)
