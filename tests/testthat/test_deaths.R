context("deaths")

test_that("no deaths occur with zero case fatality ratio", {
  set.seed(1) ## to do: replicate with different random numbers
  # single country, single age group, single risk group, all people in S
  
  ## LIFE HISTORY PARAMETER INPUTS
  ## R_0, recovery time and latent period
  life_history_params <- list(R0=1.8, TR=2.6, LP = 1.5)
  ## vaccine efficacy and initial vaccinated proportion
  # this example roughly brings effective R to 1.2
  vax_params <- list(efficacy = 1 - 1.2/1.8, propn_vax0 = 0)
  
  ## example parameters for vaccine production (see cum_vax_pool_func_closure)
  vax_production_params <- list(dummy = 0)
  
  ## example parameters for vaccine allocation
  vax_allocation_params <- list(priorities = NULL)
  
  popn_size <- 100000
  n_ages <- 1
  age_propns <- 1
  n_countries <- 1
  n_riskgroups <- 1
  
  risk_propns <- diag(1)
  risk_factors <- 1
  
  tmp <- setup_populations(popn_size,n_countries,age_propns,
                           risk_propns, risk_factors)
  
  X <- tmp$X
  labels <- tmp$labels
  
  K <- diag(1)
  latitudes <- 0
  
  C3 <- diag(1)
  
  user_specified_vax_alloc_func <- function(sum_age_risk_func, 
                                            travel_matrix,
                                            vax_allocation_params,
                                            S, E, I, R, vax_pool) {
    return(0)
  }
  
  user_specified_cum_vax_pool_func <- function(vax_production_params, t) {
    return(0)
  }
  
  ## process the vaccine production function
  cum_vax_pool_func <- cum_vax_pool_func_closure(user_specified_cum_vax_pool_func, vax_production_params)
  
  ## process the vaccine allocation function
  
  vax_allocation_func <- vaccine_allocation_closure(user_specified_vax_alloc_func, K, vax_allocation_params, labels)
  
  seed_vec <- 10

  ## gather simulation parameters
  sim_params <- list(n_countries=n_countries,
                     n_ages=n_ages,
                     n_riskgroups=n_riskgroups,
                     seed_vec = seed_vec)
  
  simulation_flags <- list(ageMixing=TRUE,
                           riskGroups=TRUE,
                           normaliseTravel=TRUE,
                           spatialCoupling=TRUE,
                           real_data = FALSE,
                           country_specific_contact = FALSE,
                           seasonal = FALSE,
                           rng_seed = 0)
  
  ## run simulation for tmax days
  tmax <- 100
  ## tdiv timesteps per day
  tdiv <- 24 ##DH

  vax_alloc_period <- 24 * 7
  
  case_fatality_ratio_vec <- 0
  
  set.seed(1)
  ## run simulation
  res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                        case_fatality_ratio_vec, X, labels, C3, K, latitudes, 
                        cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period)
  
  get_deaths_recovered <- function(res) {
    tend <- ncol(res$S)
    deaths <- as.numeric(X - res$S[,tend] - res$SV[,tend] - res$E[,tend] - res$EV[,tend] -
      res$I[,tend] - res$IV[,tend] - res$R[,tend] - res$RV[,tend])
    recovered <- unname(res$R[,tend] + res$RV[,tend])
    list("deaths" = deaths, "recovered" = recovered)
  }

  deaths_recovered <- get_deaths_recovered(res)
  expect_equal(sum(deaths_recovered$deaths), 0)
  
  case_fatality_ratio_vec <- 1
  
  set.seed(1)
  ## run simulation
  res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                        case_fatality_ratio_vec, X, labels, C3, K, latitudes, 
                        cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period)
  
  deaths_recovered2 <- get_deaths_recovered(res)
  expect_equal(sum(deaths_recovered2$recovered), 0)
  expect_equal(deaths_recovered$recovered, deaths_recovered2$deaths)

})