library(reshape2)
library(ggplot2)
library(Matrix)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=4)

# uncomment lines at bottom to run

wd <- "~/Documents/vaxedemic/" 
devtools::load_all(wd)
setwd(wd)
run_simulation_vary_vaccine <- function(efficacy, production_delay, production_rate, user_specified_vax_alloc_func) {
###################################################
## RUN CONTROL
###################################################
## Number of simulations to run
n_runs <- 4
## run simulation for tmax days
tmax <- 365
## tdiv timesteps per day
tdiv <- 24
## Seasonality resolution
seasonality_resolution <- tmax*tdiv/12 # Average seasonality into 12 blocks of time
###################################################


###################################################
## LIFE HISTORY PARAMETER INPUTS
###################################################
## R_0, recovery time and latent period
life_history_params <- list(R0=1.4, TR=2.6, LP = 1.5, case_fatality_ratio = rep(2e-2,2))

## travel parameters: scaling of off-diagonals
travel_params <- list(epsilon = 1e-3)
###################################################


###################################################
## SIMULATION OPTIONS
###################################################
simulation_flags <- list(ageMixing=TRUE,
                         riskGroups=TRUE,
                         normaliseTravel=TRUE,
                         spatialCoupling=TRUE,
                         real_data = TRUE,
                         country_specific_contact = TRUE,
                         seasonal = TRUE,
                         rng_seed = 0)
###################################################


###################################################
## VACCINE INPUTS
###################################################
## vaccine efficacy and initial vaccinated proportion

vax_params <- list(efficacy = efficacy, propn_vax0 = 0)

## vaccine production function
user_specified_cum_vax_pool_func <- produce_vax_linear_with_delay
## parameters for vaccine production (see cum_vax_pool_func_closure)
vax_production_params <- list(detection_delay = 0, production_delay = production_delay, 
                              production_rate = production_rate, max_vax = Inf)

## parameters for vaccine allocation
vax_allocation_params <- list(priorities = NULL, coverage = NULL)

## allocate and distribute vaccine every vac_alloc_period time divisions
## i.e. in this example, every 7 days
vax_alloc_period <- 24 * 7
###################################################


###################################################
## SETUP SEEDING
###################################################
# set random number generation seed
if(!is.null(simulation_flags[["rng_seed"]])) {
  set.seed(simulation_flags[["rng_seed"]])
}

## Seeding setting
if(simulation_flags[["real_data"]]) {
  seedCountries <- "China"
} else {
  seedCountries <- 1
}

## number of exposed individuals in each seeded country
seedSizes <- c(20)
## which age group(s) to seed
seedAges <- 3
## which risk group(s) to seed
seedRiskGroups <- 1


###################################################
## SETUP INPUTS
###################################################
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

#construct vector of number of exposed individuals in each location, age, risk group 
## for now, can only seed in one location, age, risk group. Vectorise later
seed_vec <- double(length(popns))
seed_vec[(which(labels$Location == seedCountries & labels$Age == seedAges &
              labels$RiskGroup == seedRiskGroups))[1]] <- seedSizes

## process the vaccine production function
cum_vax_pool_func <- cum_vax_pool_func_closure(user_specified_cum_vax_pool_func, vax_production_params)

if ("coverage" %in% names(vax_allocation_params)) {
  if (simulation_flags[["real_data"]]) {
    coverage_filename <- "data/coverage_data_intersect.csv"
    coverage <- read_coverage_data(coverage_filename, labels)
  } else {
    # every country has same seasonal coverage
    # the below line is incorrect: ada to fix
    # coverage <- rep(1/n_countries, n_countries)#
    stop("uniform coverage not yet implemented")
  }
  vax_allocation_params$coverage <- coverage
}

## process the vaccine allocation function
vax_allocation_func <- vaccine_allocation_closure(user_specified_vax_alloc_func,
                                                  travelMatrix, vax_allocation_params, labels)

## gather simulation parameters
sim_params <- list(n_countries=n_countries,
                   n_ages=n_ages,
                   n_riskgroups=n_riskgroups,
                   seed_vec = seed_vec,
                   seasonality_resolution=seasonality_resolution,
                   tdelay=180,
                   amp=.7)
###################################################
res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                      case_fatality_ratio_vec, popns, labels, contactMatrix, travelMatrix, latitudes, 
                      cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period,
                      n_runs=n_runs)

summary_stats <- deaths_GAR_df(res, popns)
mean_summary_stats <- apply(summary_stats, 2, mean)
country_attack_rate <- lapply(res, country_attack, popns, labels)
list(country_attack_rate = country_attack_rate, 
     labels = labels, 
     summary_stats = mean_summary_stats)
}

run_many_simulations_vary_vaccine <- function(vax_strategy) {
  production_delay <- c(0, 365/4, 365/2)
  pars <- expand.grid(efficacy = c(0.7,1), 
                      production_delay = production_delay, 
                      production_rate = c(550e6/(365/4), 7e9))
  pars <- rbind(data.frame(efficacy = 0, production_delay = 0, production_rate = 0), pars)

  if(vax_strategy == "incidence") {
    user_specified_vax_alloc_func <- vaccinate_by_incidence
  } else if (vax_strategy == "current_seasonal_alloc") {
    user_specified_vax_alloc_func <- vaccinate_by_current_seasonal_alloc
  } else {
    stop("unknown vaccination strategy")
  }
  
  results <- lapply(seq_len(nrow(pars)), function(x) run_simulation_vary_vaccine(pars[x, 1], pars[x, 2], pars[x, 3],
                                                                                 user_specified_vax_alloc_func))

  labels <- results[[1]]$labels
  country_attack_rate <- lapply(results, function(x) x$country_attack_rate)
  g <- lapply(country_attack_rate, plot_country_attack_rates, labels)
  Map(function(x, y) ggsave(paste0("country_attack_vaccinate_by_", vax_strategy, x,".pdf"), y$plot,
                            width = 21, height = 29, units = "cm"), seq_along(g), g)


  global_attack <- vapply(results, function(x) x$summary_stats[["global_attack"]], double(1))
  summary_stats <- cbind(pars, data.frame(global_attack = global_attack))
  saveRDS(summary_stats, paste0("summary_stats_vax_by_", vax_strategy, ".rds"))
  summary_stats$production_delay <- as.factor(summary_stats$production_delay)
  summary_stats$production_rate <- as.factor(summary_stats$production_rate)
  levels(summary_stats$production_rate) <- c("0", "current", "large")
  levels(summary_stats$production_delay) <- c(0, 3, 6)
  g <- ggplot(summary_stats[-1,]) +
    geom_bar(aes(x = production_delay, y = global_attack,
                 group = production_rate, fill = production_rate),
             stat = "identity", position = "dodge") +
    facet_grid(.~efficacy) +
    xlab("Production delay (months)") +
    ylab("Global attack rate") +
    geom_hline(yintercept = summary_stats[1,"global_attack"], lty = "dotted")
  ggsave(paste0("vaccinate_by_", vax_strategy, ".pdf"), g)
  g
}

# uncomment to run
# g <- run_many_simulations_vary_vaccine("incidence")
# g <- run_many_simulations_vary_vaccine("current_seasonal_alloc")
