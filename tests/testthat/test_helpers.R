context("helpers")

# test_that("round_preserve_sum works", {
#   n <- 10
#   vec <- runif(n, min = 0, max = 5)
#   vec[1] <- ceiling(sum(vec[-1])) - sum(vec[-1])
#   sum_vec <- sum(vec)
#   expect_equal(sum_vec, round(sum_vec))
#   n_reps <- 1e5
#   replicate_round_preserve_sum <- replicate(n_reps, round_preserve_sum(vec))
#   expect_true(all(colSums(replicate_round_preserve_sum) == sum_vec))
#   n_rounded_up <- rowSums(replicate_round_preserve_sum) - floor(vec) * n_reps
#   probs <- vec - floor(vec)
#   expect_equal(order(probs), order(n_rounded_up))
# })

test_that("gen_int_sum_int works", {
  set.seed(1)
  n_trials <- 100
  n_int <- floor(runif(n_trials, min = 1, max = 10))
  sum_int <- floor(runif(n_trials, min = 1, max = 100))
  gen_int <- Map(gen_int_sum_int, n_int, sum_int)
  
  expect_equal(vapply(gen_int, length, double(1)), n_int)
  expect_equal(vapply(gen_int, sum, double(1)), sum_int)
})

test_that("distribute_vax_among_age_risk_closure works", {
  set.seed(1) ## to do: replicate with different random numbers
  # single country, single age group, single risk group, all people in S
  
  popn_size <- 100000
  n_ages <- 1
  age_propns <- 1
  n_countries <- 1
  n_riskgroups <- 1
  
  risk_propns <- rep(1/n_riskgroups,n_riskgroups) ## Assume risk groups are uniformly distributed
  risk_propns <- matrix(rep(risk_propns,each=n_ages),ncol=n_riskgroups) ## Assume that proportion of ages in each risk group are the same for all ages
  risk_factors <- rep(1, n_riskgroups) ## Assume that each risk group has same modifier
  
  tmp <- setup_populations(popn_size,n_countries,age_propns,
                           risk_propns, risk_factors)
  
  X <- tmp$X
  labels <- tmp$labels
  
  sum_age_risk <- sum_age_risk_closure(labels)
  
  n_vax_allocated <- floor(runif(1, max = popn_size))
  
  n_states <- 4
  S <- popn_size
  E <- I <- R <- 0

  priorities_base <- expand.grid("RiskGroup" = seq_len(n_riskgroups), "Age" = seq_len(n_ages))
  priorities <- NULL
  distribute_vax_among_age_risk <- distribute_vax_among_age_risk_closure(priorities, labels)
  vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
  
  # expect everything to be allocated to S
  expect_equal(as.numeric(vax_alloc), c(n_vax_allocated,0,0,0))
  
  # single country, single age group, single risk group, 
  # people distributed among S, E, I, R
  pop_sizes_SEIR <- gen_int_sum_int(n_states, popn_size)
  
  S <- pop_sizes_SEIR[1]
  E <- pop_sizes_SEIR[2]
  I <- pop_sizes_SEIR[3]
  R <- pop_sizes_SEIR[4]
  
  n_vax_allocated <- popn_size
  vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
  
  # expect everyone to get vaccinated in each infection state
  expect_equal(as.numeric(vax_alloc), c(S,E,I,R))
  
  n_vax_allocated <- floor(runif(1, max = popn_size))
  vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
  
  # expect proportional distribution of vaccines
  vax_alloc_vec <- as.numeric(vax_alloc)
  fractional_expected_vax_alloc <- n_vax_allocated * c(S,E,I,R) / popn_size
  expect_equal(sum(vax_alloc_vec), n_vax_allocated)
  expect_true(all(vax_alloc_vec >= floor(fractional_expected_vax_alloc)))
  expect_true(all(vax_alloc_vec <= ceiling(fractional_expected_vax_alloc)))
  
  # more than one country, single age group, single risk group, 
  # n_vax_allocated in each country, people all in S
  
  popn_size <- c(1e5, 2e5)
  n_countries <- length(popn_size)
  
  X <- popn_size
  labels <- cbind(X,expand.grid("Location"=1:n_countries, "RiskGroup"=1:n_riskgroups,"Age"=1:n_ages))

  sum_age_risk <- sum_age_risk_closure(labels)
  
  n_vax_allocated <- floor(runif(n_countries, min = 0, max = popn_size))
  
  S <- popn_size
  E <- I <- R <- double(n_countries)

  distribute_vax_among_age_risk <- distribute_vax_among_age_risk_closure(priorities, labels)
  vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)

  # expect all vaccines to be distributed to S in each country
  expect_equal(unname(unlist(vax_alloc)), c(n_vax_allocated,E,I,R))
  
  # single country, more than one age/risk group, people all in S, priority is NULL
  
  popn_size <- 1e5
  n_countries <- length(popn_size)
  
  n_ages <- 3
  age_propns <- normalise(runif(n_ages))
  n_countries <- 1
  n_riskgroups <- 2
  
  risk_propns <- rep(1/n_riskgroups,n_riskgroups) ## Assume risk groups are uniformly distributed
  risk_propns <- matrix(rep(risk_propns,each=n_ages),ncol=n_riskgroups) ## Assume that proportion of ages in each risk group are the same for all ages
  risk_factors <- rep(1, n_riskgroups) ## Assume that each risk group has same modifier
  
  tmp <- setup_populations(popn_size,n_countries,age_propns,
                           risk_propns, risk_factors)
  
  X <- tmp$X
  labels <- tmp$labels
  
  n_vax_allocated <- popn_size

  S <- X
  E <- I <- R <- S * 0

  distribute_vax_among_age_risk <- distribute_vax_among_age_risk_closure(priorities, labels)
  vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
  
  # expect n_vax_allocated spread out proportionally
  expect_equal(unname(unlist(vax_alloc)), c(S,E,I,R))
  
  n_vax_allocated <- floor(runif(1, max = popn_size))
  vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
  # expect all vaccines to end up somewhere
  expect_equal(sum(unname(unlist(vax_alloc))), n_vax_allocated)
  # expect n_vax_allocated spread out proportionally
  fractional_expected_vax_alloc <- n_vax_allocated * S / popn_size
  expect_true(all(vax_alloc$S >= floor(fractional_expected_vax_alloc)))
  expect_true(all(vax_alloc$S <= ceiling(fractional_expected_vax_alloc)))
})


  # # single country, more than one age/risk group, people all in S, priority is specified
  # 
  # priorities_base <- expand.grid("RiskGroup"=seq_len(n_riskgroups),"Age"=seq_len(n_ages))
  # priorities <- cbind(priorities_base, data.frame("Priority" = seq_len(n_riskgroups * n_ages)))
  # 
  # distribute_vax_among_age_risk <- distribute_vax_among_age_risk_closure(priorities, labels)
  # 
  # # first test naive allocation where number of vaccines allocated in age/risk 
  # # group can exceed the number of individuals in that age/risk group
  # fractional_expected_vax_alloc <- n_vax_allocated * normalise(X * priorities$Priority)
  # vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
  # 
  # expect_equal(sum(unname(unlist(vax_alloc))), n_vax_allocated)
  # expect_true(all(vax_alloc$S >= floor(fractional_expected_vax_alloc)))
  # expect_true(all(vax_alloc$S <= ceiling(fractional_expected_vax_alloc)))
  # 
  # priorities <- cbind(priorities_base, data.frame("Priority" = rep(1e-5, n_riskgroups * n_ages)))
  # priorities[1,"Priority"] <- 1
  # 
  # distribute_vax_among_age_risk <- distribute_vax_among_age_risk_closure(priorities, labels)
  # fractional_expected_vax_alloc <- n_vax_allocated * normalise(X * priorities$Priority)
  # vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)

  # # while(!isTRUE(all.equal(vax_alloc, actual_alloc))) {
  # 
  #   # vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
  #   # expected_vax_alloc <- vax_pool * normalise(X * priorities$Priority)
  # #   actual_alloc <- Map(pmin, vax_alloc, list(S, E, I, R))
  # #   browser()
  # #   expected_actual_alloc <- Map(pmin, expected_vax_alloc, list(S, E, I, R))
  # #   sum_expected_alloc <- Map(`+`, sum_expected_alloc, expected_actual_alloc)
  # # 
  # #   ## update current vax pool
  # #   vax_pool <- vax_pool - sum(unlist(actual_alloc))
  # #   
  # #   sum_vax_alloc <- Map(`+`, sum_vax_alloc, actual_alloc)
  # # 
  # #   S <- S - actual_alloc$S
  # #   E <- E - actual_alloc$E
  # #   I <- I - actual_alloc$I
  # #   R <- R - actual_alloc$R
  # # }
  # 
  # # expect all vaccines to end up somewhere
  # expect_equal(sum(unname(unlist(sum_vax_alloc))), n_vax_allocated)
  # fractional_expected_vax_alloc <- n_vax_allocated * normalise(X * priorities$Priority)
  # browser()
  # # expect n_vax_allocated spread out proportionally
  # expect_true(all(sum_vax_alloc$S >= floor(fractional_expected_vax_alloc)))
  # expect_true(all(sum_vax_alloc$S <= ceiling(fractional_expected_vax_alloc)))
  # # expect n_vax_allocated spread out proportionally
  
