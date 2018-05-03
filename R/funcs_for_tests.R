#' function to simulate the final size distribution for a single country
#' 
#' function to simulate thefinal size distribution for a single country, for a 
#' give effective reproduction number.
#' 
#' @param R_eff numeric vector of length 1: effective reproduction number.
#' R_0 = 1.4
#' @param change_efficacy logical: if TRUE, achieve R_eff by vaccinating all
#' individuals but having vaccine efficacy less than unity.  If FALSE, achieve
#' R_eff by vaccinating a proportion of individual, with vaccine efficacy equal
#' to unity.
#' @return number vector of length 100: final sizes of 100 simulations.
get_final_sizes_with_R_eff <- function(R_eff, change_efficacy) {
  
  ###################################################
  ## RUN CONTROL
  ###################################################
  ## Number of simulations to run
  n_runs <- 100
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
  life_history_params <- list(R0=1.4, TR=2.6, LP = 1.5, case_fatality_ratio = 0)
  
  ## travel parameters: scaling of off-diagonals
  travel_params <- list(epsilon = 1e-3)
  ###################################################
  
  
  ###################################################
  ## SIMULATION OPTIONS
  ###################################################
  simulation_flags <- list(ageMixing=FALSE,
                           riskGroups=FALSE,
                           normaliseTravel=TRUE,
                           spatialCoupling=FALSE,
                           real_data = FALSE,
                           country_specific_contact = FALSE,
                           seasonal = FALSE,
                           rng_seed = 0)
  ###################################################
  
  
  ###################################################
  ## VACCINE INPUTS
  ###################################################
  ## vaccine efficacy and initial vaccinated proportion
  
  # achieve R_eff either by having all individuals vaccinated at start, with vaccine
  # efficacy below 1
  if(change_efficacy) {
    vax_params <- list(efficacy = 1 - R_eff/life_history_params$R0, propn_vax0 = 1)
  } else {
    # or have a proportion of individuals vaccinated at start, with vaccine
    # efficacy = 1
    vax_params <- list(efficacy = 1, propn_vax0 = R_eff/life_history_params$R0)
  }
  
  
  ## vaccine production function
  user_specified_cum_vax_pool_func <- produce_vax_linear_with_delay
  ## parameters for vaccine production (see cum_vax_pool_func_closure)
  vax_production_params <- list(detection_delay = 0, production_delay = 0, 
                                production_rate = 0, max_vax = 0)
  
  ## vaccination strategy
  user_specified_vax_alloc_func <- vaccinate_by_incidence
  ## parameters for vaccine allocation
  vax_allocation_params <- list(priorities = NULL)
  
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
  seedCountries <- 1
  
  ## number of exposed individuals in each seeded country
  seedSizes <- c(20)
  ## which age group(s) to seed
  seedAges <- 1
  ## which risk group(s) to seed
  seedRiskGroups <- 1
  
  
  ###################################################
  ## SETUP INPUTS
  ###################################################
  all_inputs <- setup_sim_data(simulation_flags,
                               life_params,
                               travel_params,
                               popn_size=100000,
                               n_riskgroups=1,
                               n_ages=1,
                               age_propns=1,
                               n_countries=1)
  
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
  
  final_sizes <- summary_stats$global_attack * popns[1,1]
}

#' Britton2010 10.1016/j.mbs.2010.01.006
calc_approx_final_size_distribution <- function(R_0, pop_size) {
  z_star <- calc_final_size_propn_deterministic(R_0)
  mu <- z_star * pop_size
  sigma <- sqrt(pop_size * z_star * (1 - z_star) * (1 + (1 - z_star) * R_0^2)) / 
    (1 - (1 - z_star) * R_0)
  list("mu" = mu, "sigma" = sigma)
}

calc_final_size_propn_deterministic <- function(R_0, S_0 = 1) {
  final_size_equation <- function(z) {
    1 - z - S_0*exp(-R_0 * z)
  }
  
  epsilon <- 1e-10
  uniroot(final_size_equation, c(epsilon, 1 - epsilon))$root
}

two_tailed_p_value <- function(q, mu, sigma) {
  p <- pnorm(q, mu, sigma)
  p <- min(p, 1 - p)
  2 * p
}

z_test <- function(q, mu, sigma) {
  p_values <- vapply(q, function(q) two_tailed_p_value(q, mu, sigma), double(1))
  n_tests <- length(q)
  list("min_p" = min(p_values), "n_tests" = n_tests)
} 

print_z_test_result <- function(z_test_result, n_tests = 1) {
  format_results <- function(x) {
    formatC(x, format = "e", digits = 2)
  }
  print(paste0("z test p-value = ", format_results(z_test_result$min_p)))
  print(paste0("z test n_tests = ", z_test_result$n_tests * n_tests))
}

result_z_test_adjust_alpha <- function(result_z_test, alpha) {
  adjusted_alpha <- alpha/result_z_test$n_tests
  significant <- result_z_test$min_p < adjusted_alpha
  significant
}
