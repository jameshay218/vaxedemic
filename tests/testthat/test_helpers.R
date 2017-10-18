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

test_that("distribute_vax_among_age_risk_closure works", {
  
  # set up example
  
  popn_size <- 100000
  ## Setup age propns
  n_ages <- 3
  age_propns <- c(5,14,45+16)/80
  n_countries <- 10
  n_riskgroups <- 2
  n_groups <- n_ages * n_riskgroups
  
  risk_propns <- rep(1/n_riskgroups,n_riskgroups) ## Assume risk groups are uniformly distributed
  risk_propns <- matrix(rep(risk_propns,each=n_ages),ncol=n_riskgroups) ## Assume that proportion of ages in each risk group are the same for all ages
  risk_factors <- rep(1, n_riskgroups) ## Assume that each risk group has same modifier
  
  ## construct demography matrix
  tmp <- setup_populations(popn_size,n_countries,age_propns, 
                           risk_propns, risk_factors)
  
  X <- tmp$X
  labels <- tmp$labels
  

  
  sum_age_risk <- sum_age_risk_closure(labels)
  
  vax_alloc_scaling <- 1
  n_vax_allocated <- round(runif(n_countries) * sum_age_risk(X) * vax_alloc_scaling)
  
  n_states <- 4
  SEIR <- matrix(runif(n_countries * n_groups * n_states), ncol = n_states)
  ## make SEIR in each location, country, age group sum to population size
  SEIR <- sweep(SEIR, 1, rowSums(SEIR) / X, FUN = "/")
  SEIR <- round(SEIR)
  
  S <- SEIR[,1]
  E <- SEIR[,2]
  I <- SEIR[,3]
  # ensure SEIR sum to X after rounding
  R <- as.numeric(X - S - E - I)

  priorities_base <- expand.grid("RiskGroup" = seq_len(n_riskgroups), "Age" = seq_len(n_ages))
  priorities <- cbind(priorities_base, data.frame("Priority" = rep(1 / n_groups, n_groups)))
  distribute_vax_among_age_risk <- distribute_vax_among_age_risk_closure(priorities, labels)
  vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)

  expect_equal(sum_age_risk(vax_alloc$S + vax_alloc$E + vax_alloc$I + vax_alloc$R), n_vax_allocated)
  
  priorities <- cbind(priorities_base, data.frame("Priority" = seq_len(n_groups)))

  distribute_vax_among_age_risk <- distribute_vax_among_age_risk_closure(priorities, labels)
  vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
  expect_equal(sum_age_risk(vax_alloc$S + vax_alloc$E + vax_alloc$I + vax_alloc$R), n_vax_allocated)
  
  priorities <- cbind(priorities_base, data.frame("Priority" = c(1, double(n_groups - 1))))
  
  distribute_vax_among_age_risk <- distribute_vax_among_age_risk_closure(priorities, labels)
  vax_alloc <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
  expect_equal(sum_age_risk(vax_alloc$S + vax_alloc$E + vax_alloc$I + vax_alloc$R), n_vax_allocated)
  
})