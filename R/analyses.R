#' example user-specified vaccine allocation function
#' 
#' allocate vaccines according to absolute incidence in each location, then
#' allocate uniformly within each country (without discriminating between
#' ages/risk groups/infection statuses)
#' needs to take the arguments
#' @param sum_age_risk_func a function which sums a vector across age and risk groups (created automatically in the code)
#' @param travel_matrix a square matrix with side length n_countries. Travel_matrix[x,y] is the proportion of time an individual in location x spends in location y
#' @param vax_allocation_params list of parameters for vaccine allocation
#' @param S numeric vector of length n. number of unvaccinated susceptibles.
#' @param n is number of locations * number of age groups * number of risk groups
#' @param E numeric vector of length n. number of unvaccinated exposed
#' @param I numeric vector of length n. number of unvaccinated infectious
#' @param R numeric vector of length n. number of unvaccinated recovered
#' @param vax_pool numeric vector of length 1. number of vaccines available needs to return
#' @param n_vax_allocated numeric vector of length n_countries.  number of vaccines allocated to each country
#'
#' @export
user_specified_vax_alloc_func <- function(sum_age_risk_func, 
                                          travel_matrix,
                                          vax_allocation_params,
                                          S, E, I, R, vax_pool) {
  # find incidence in each country
  # incidence proportional to E
  incidence_by_country <- sum_age_risk_func(E)
  if(any(incidence_by_country > 0)) {
    n_vax_allocated <- incidence_by_country / sum(incidence_by_country) * vax_pool
  } else {
    n_vax_allocated <- incidence_by_country * 0 # allocate nothing
  }
  return(n_vax_allocated)
}

#' example function of vaccine production:
#' 
#' no vaccine produced until time vax_production_params[["detection_delay"]] + vax_production_params[["production_delay"]], then constant production rate until max number of doses ever made reached, then no production
#' @param vax_production_params list of parameters for vaccine production
#' @param t scalar time
#' @return needs to return a scalar: the number of vaccines ever produced up to time t
#' @export
user_specified_cum_vax_pool_func <- function(vax_production_params, t) {
  t_since_production <- t - (vax_production_params[["detection_delay"]] + 
                               vax_production_params[["production_delay"]])
  if(t_since_production < 0) {
    0
  } else {
    min(vax_production_params[["max_vax"]],
        t_since_production * vax_production_params[["production_rate"]])
  }
}
