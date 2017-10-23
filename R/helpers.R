#' returns a function to sum a state vector (e.g. the current number of susceptibles in each
#' country, age and risk group) across age and risk groups
#' 
#' @param labels data frame with n rows, and columns "Location", "RiskGroup" and "Age"
#' @return a function to sum a state vector (e.g. the current number of susceptibles in each
#' country, age and risk group) across age and risk groups
#' @export
sum_age_risk_closure <- function(labels) {
  function(state) {
    aggregate(state, by = list(labels$Location), FUN = sum)[,2]
  }
}

#' returns a function which, 
#' given the number of vaccines allocated to each country, distributes vaccines
#'  between age and risk groups
#'  
#'  @param priorities data frame with rows n_groups = n_ages * n_riskgroups.
#'  Contains the columns "RiskGroup","Age" and "Priority".
#'  If x vaccines are allocated to a country, and the unvaccinated population in
#'  each age/risk group in that country is given by a vector p of length n_groups,
#'  then the number of vaccines allocated to each age/risk group is
#'  x * p * y / sum(p * y) where y is a normalised version of priorities$Priority.
#'  For example, to give vaccines only to risk/age group 1, y = [1,0,0,...].
#'  For non-discriminatory distribution, y = [1,1,1,...] / n_groups.
#'  Note that if y_i = 0 for any age/risk group i, vaccines will never be allocated
#'  to that age/risk group even if all other age/risk groups have been fully vaccinated --
#'  best to use y_i = epsilon << 1 instead.
#'  For consistency with labels, 
#'  @param labels data frame with n rows, and columns "Location", "RiskGroup" and "Age".
#'  n = n_countries * n_ages * n_riskgroups.
#'  @return a function f.
#'  f takes the arguments
#'  n_vax_allocated: numeric vector of length n_countries. 
#'  The number of vaccines allocated to each country.
#'  n_unvaccinated: numeric vector of length n.
#'  The number of unvaccinated individuals in each location, age, risk group.
#'  f returns a numeric vector of length n: the
#'  number of vaccines allocated to each country, age, risk group
distribute_vax_among_age_risk_closure <- function(priorities, labels) {
  
  ## if priorities not given, default behaviour is non-discriminatory distribution
  if(is.null(priorities)) {
    n_groups <- length(unique(labels$Age)) * length(unique(labels$RiskGroup))
    n_countries <- round(nrow(labels) / n_groups)
    y_long_vec <- rep(1 / n_groups, nrow(labels))
  } else {
    stop("code with specified priorities doesn't work yet")
    ## sort priorities
    # priorities <- priorities[order(priorities$Age, priorities$RiskGroup),]
    # n_groups <- nrow(priorities)
    # n_countries <- round(nrow(labels) / n_groups)
    # y <- normalise(priorities$Priority)
    # y_long_vec <- rep(y, each = n_countries)
  }
  
  state_names <- c("S", "E", "I", "R")
  n_states <- length(state_names)
  ## non-discriminatory distribution among S, E, I, R
  y_long_vec <- rep(y_long_vec, times = n_states)

  sum_age_risk <- sum_age_risk_closure(labels)
  
  # index of each location in y_long_vec
  loc_idx <- lapply(seq_len(n_countries), function(x) seq(x, length(y_long_vec), by = n_countries))
  
  apply_by_loc_closure <- function(loc_idx, n_countries) {
    f <- function(f, vec) {
      list_by_loc <- lapply(loc_idx, function(x) f(vec[x]))
      return(as.numeric(matrix(unlist(list_by_loc), nrow = n_countries, byrow = TRUE)))
    }
    f
  }
  
  apply_by_loc <- apply_by_loc_closure(loc_idx, n_countries)
  
  f <- function(n_vax_allocated, S, E, I, R) {
    
    # ensure that an integer number of vaccines is allocated to each country
    n_vax_allocated <- round_preserve_sum(n_vax_allocated)
    
    # distribute vaccines between age and risk groups
    vax_alloc <- rep(n_vax_allocated, times = n_groups * n_states)
    distribution_factor <- c(S, E, I, R) * y_long_vec
    distribution_factor <- apply_by_loc(normalise, distribution_factor)
    vax_alloc <- vax_alloc * distribution_factor

    # round ensuring that the number of vaccines allocated to each country stays the same  
    vax_alloc <- apply_by_loc(round_preserve_sum, vax_alloc)

    # split the vaccine allocation vector into S, E, I, R
    vax_alloc <- split(vax_alloc, ceiling(seq_along(vax_alloc) / n_groups / n_countries))
    names(vax_alloc) <- state_names
    return(vax_alloc)
  }
  f
}

#' rounds a numeric vector probabilistically such that its sum (an integer) is preserved
#' 
#' @param vec numeric vector which sums to an integer
#' @return numeric vector of the same length as vec: rounded vector
round_preserve_sum <- function(vec) {
  
  # if vector already integers, do nothing
  if(isTRUE(all.equal(vec, round(vec)))) {
    return(vec)
  }
  
  sum_vec <- sum(vec)
  rounded_sum <- round(sum_vec)
  if(!isTRUE(all.equal(rounded_sum, sum_vec))) {
    stop("vector passed to round_preserve_sum does not sum to integer")
  }
  
  n_round_up <- rounded_sum - sum(floor(vec))
  round_up_idx <- sample.int(length(vec), size = n_round_up, replace = FALSE,
                             prob = vec - floor(vec))
  vec_out <- floor(vec)
  vec_out[round_up_idx] <- ceiling(vec[round_up_idx])
  vec_out
}

#' normalises a vector so that it sums to 1
#' 
#' @param vec numeric vector to be normalised
#' @return normalised vector of same length as vec
normalise <- function(vec) {
  
  if(all(vec == 0)) {
    return(vec)
  } else if(any(vec < 0)) {
    stop("vector to be normalised must have all non-negative elements")
  }
  
  vec / sum(vec)
  
}

#' generate n_int random integers that sum to sum_int
#' 
#' @param n_int number of random integers to generate
#' @param sum_int integer to which those integers should sum to
#' @return numeric vector of length n_int
gen_int_sum_int <- function(n_int, sum_int) {
  round_preserve_sum(normalise(runif(n_int)) * sum_int)
}

# a closure to make a function which returns a vector of length
# n = n_countries * n_ages * n_risk_groups:
# the number of vaccines to be allocated to each location, age, risk group
#'
#' @param user_specified_vax_alloc_func function to allocate vaccines
#' given the arguments
#' #' @param sum_age_risk_func: a function which sums a vector across age and risk groups
#' #' (created in the closure)
#' #' @param travel_matrix: see below
#' #' @param vax_allocation_params: see below
#' #' @param S: numeric vector of length n. number of unvaccinated susceptibles
#' #' @param E: numeric vector of length n. number of unvaccinated exposed
#' #' @param I: numeric vector of length n. number of unvaccinated infectious
#' #' @param R: numeric vector of length n. number of unvaccinated recovered
#' #' @param vax_pool: numeric vector of length 1. number of vaccines available
#' @param travel_matrix square matrix with side length n_countries. 
#' travel_matrix[x,y] is the proportion of time an individual in location x spends
#' in location y
#' @param vax_allocation_params list of parameters for vaccine allocation
#' @param labels data frame containing location, age and risk group corresponding to a vector element
#' @return a function which returns a vector of length
# n_countries * n_ages * n_risk_groups:
# the number of vaccines to be allocated to each location, age, risk group
vaccine_allocation_closure <- function(user_specified_vax_alloc_func, 
                                       travel_matrix, vax_allocation_params, labels) {
  
  ## create a function which sums a state vector across age and risk groups
  sum_age_risk_func <- sum_age_risk_closure(labels)
  distribute_vax_among_age_risk <- 
    distribute_vax_among_age_risk_closure(vax_allocation_params$priorities, labels)
  
  function(S, E, I, R, vax_pool) {
    
    # vax_pool can be not an integer -- round down
    vax_pool <- floor(vax_pool)
    
    # allocate vaccines to countries
    n_vax_allocated <- user_specified_vax_alloc_func(sum_age_risk_func, travel_matrix,
                                                     vax_allocation_params,
                                                     S, E, I, R, vax_pool)
    
    # distribute vaccines in each country between individuals of different
    # ages, risk groups, and infection statuses
    split_allocation <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
    return(split_allocation)
  }
}

#' a closure to make a function which returns a scalar:
# the number of vaccines which have ever existed at time t
# (may involve integrating the production rate curve if you want,
# otherwise provide directly)
#' @param user_specified_cum_vax_pool_func: function of vax_production_params and t
#' which returns the number of vaccines which have ever existed at time t
#' @param vax_production_params: list of parameters for vaccine production
#' @return function of t which returns
#' the number of vaccines which have ever existed at time t
cum_vax_pool_func_closure <- function(user_specified_cum_vax_pool_func, vax_production_params) {
  f <- function(t) {
    user_specified_cum_vax_pool_func(vax_production_params, t)
  }
  f
}