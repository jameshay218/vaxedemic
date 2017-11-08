context("deaths")

test_that("one risk group: no deaths occur with zero case fatality ratio; 
          the number of deaths with a CFR of 1 is equal to 
          the number of recovereds with a CFR of 0", {

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
  
  case_fatality_ratio_vec <- 0
  
  approx_final_size_distribution <- calc_approx_final_size_distribution(life_history_params$R0, popn_size)
  
  set.seed(1)
  ## run simulation
  ## run simulation for tmax days
  ## tdiv timesteps per day
  args <- list("simulation_flags" = simulation_flags, 
               "life_history_params" = life_history_params, 
               "vax_params" = vax_params, 
               "sim_params" = sim_params,
  "case_fatality_ratio_vec" = case_fatality_ratio_vec, 
  "X" = X, 
  "labels" = labels, 
  "C3" = C3, 
  "K" = K, 
  "latitudes" = latitudes, 
  "cum_vax_pool_func" = cum_vax_pool_func, 
  "vax_allocation_func" = vax_allocation_func, 
  "tmax" = 100, 
  "tdiv" = 24, 
  "vax_alloc_period" = 24 * 7)
  
  deaths_recovered <- do.call(run_and_get_deaths_recovered, args)

  expect_equal(sum(deaths_recovered$deaths), 0)
  
  args$case_fatality_ratio_vec <- 1
  
  set.seed(1)
  ## run simulation
  deaths_recovered2 <- do.call(run_and_get_deaths_recovered, args)
  
  expect_equal(sum(deaths_recovered2$recovered), 0)
  expect_equal(deaths_recovered$recovered, deaths_recovered2$deaths)
  
  args$case_fatality_ratio_vec <- .5
  
  n_replicates <- 10
  deaths_recovered3 <- parallel_wrapper(seq_len(n_replicates), function(x) do.call(run_and_get_deaths_recovered, args))
  
  z_test_result <- z_test(unlist(deaths_recovered3), 
                          mu = approx_final_size_distribution$mu / 2, 
                          sigma = approx_final_size_distribution$sigma / 2)
  
  print_z_test_result(z_test_result)
})

test_that("two risk groups: no deaths occur with zero case fatality ratio; 
          the number of deaths with a CFR of 1 is equal to 
          the number of recovereds with a CFR of 0",{
  # single country, single age group, two risk groups, all people in S
  
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
  
  risk_propns <- matrix(c(.3,.7),nrow = n_countries)
  n_riskgroups <- length(risk_propns)
  risk_factors <- rep(1, n_riskgroups)
  age_specific_riskgroup_factors <- matrix(rep(risk_factors,each=n_ages),
                                           ncol=n_riskgroups)
  
  tmp <- setup_populations(popn_size,n_countries,age_propns,
                           risk_propns, risk_factors)
  
  X <- tmp$X
  labels <- tmp$labels
  
  K <- diag(n_countries)
  latitudes <- double(n_countries)
  
  C1 <- diag(n_countries)
  risk <- c(t(age_specific_riskgroup_factors))
  risk_matrix <- t(kronecker(risk,matrix(1,1,n_riskgroups*n_ages)))
  C2 <- kronecker(C1, matrix(1,n_riskgroups,n_riskgroups))
  C3 <- C2*risk_matrix
  
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
  
  seed_vec <- c(10,0)
  
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
  
  case_fatality_ratio_vec <- c(0,0)
  
  args <- list("simulation_flags" = simulation_flags, 
               "life_history_params" = life_history_params, 
               "vax_params" = vax_params, 
               "sim_params" = sim_params,
               "case_fatality_ratio_vec" = case_fatality_ratio_vec, 
               "X" = X, 
               "labels" = labels, 
               "C3" = C3, 
               "K" = K, 
               "latitudes" = latitudes, 
               "cum_vax_pool_func" = cum_vax_pool_func, 
               "vax_allocation_func" = vax_allocation_func, 
               "tmax" = 100, 
               "tdiv" = 24, 
               "vax_alloc_period" = 24 * 7)
  
  set.seed(1)
  ## run simulation
  deaths_recovered <- do.call(run_and_get_deaths_recovered, args)
  expect_equal(sum(deaths_recovered$deaths), 0)
  
  args$case_fatality_ratio_vec <- c(1,1)
  
  set.seed(1)
  ## run simulation
  deaths_recovered2 <- do.call(run_and_get_deaths_recovered, args)
  
  expect_equal(sum(deaths_recovered2$recovered), 0)
  expect_equal(deaths_recovered$recovered, deaths_recovered2$deaths)
  
  args$case_fatality_ratio_vec <- c(0,1)
  deaths_recovered3 <- do.call(run_and_get_deaths_recovered, args)
  expect_equal(deaths_recovered3$deaths[1], 0)
  expect_equal(deaths_recovered3$recovered[2], 0)
  
  args$case_fatality_ratio_vec <- c(.4,.6)
  
  n_replicates <- 10
  deaths_recovered3 <- parallel_wrapper(seq_len(n_replicates), function(x) do.call(run_and_get_deaths_recovered, args))
  
  deaths <- vapply(deaths_recovered3, function(x) x[["deaths"]], double(n_riskgroups))

  approx_final_size_distribution <- calc_approx_final_size_distribution(life_history_params$R0, popn_size)

  z_test_result <- lapply(seq_len(n_riskgroups), function(x) z_test(deaths[x, ], 
         approx_final_size_distribution$mu * args$case_fatality_ratio_vec[x] * risk_propns[x],
         approx_final_size_distribution$sigma * args$case_fatality_ratio_vec[x] * risk_propns[x]))
  
  recovered <- vapply(deaths_recovered3, function(x) x[["recovered"]], double(n_riskgroups))
  
  z_test_result2 <- lapply(seq_len(n_riskgroups), function(x) z_test(recovered[x, ], 
          approx_final_size_distribution$mu * (1 - args$case_fatality_ratio_vec[x]) * risk_propns[x],
          approx_final_size_distribution$sigma * (1 - args$case_fatality_ratio_vec[x]) * risk_propns[x]))
  
  lapply(c(z_test_result, z_test_result2), function(x) print_z_test_result(x, 4))
  
})