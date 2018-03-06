#' allocate vaccines according to absolute incidence
#' 
#' allocate vaccines according to absolute incidence in each location, then
#' allocate uniformly within each country (without discriminating between
#' ages/risk groups/infection statuses)
#' 
#' @param sum_age_risk_func a function which sums a vector across age and risk groups (created automatically in the code)
#' @param travel_matrix a square matrix with side length n_countries. Travel_matrix[x,y] is the proportion of time an individual in location x spends in location y
#' @param vax_allocation_params list of parameters for vaccine allocation
#' @param S numeric vector of length n. number of unvaccinated susceptibles.
#' @param n is number of locations * number of age groups * number of risk groups
#' @param E numeric vector of length n. number of unvaccinated exposed
#' @param I numeric vector of length n. number of unvaccinated infectious
#' @param R numeric vector of length n. number of unvaccinated recovered
#' @param vax_pool numeric vector of length 1. number of vaccines available 
#' @return n_vax_allocated numeric vector of length n_countries.  number of vaccines allocated to each country
#' @export
vaccinate_by_incidence <- function(sum_age_risk_func, 
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


#' allocate vaccines according to current seasonal allocation
#' 
#' allocate vaccines according to current seasonal allocation in each location, then
#' allocate uniformly within each country (without discriminating between
#' ages/risk groups/infection statuses)
#'
#' @param sum_age_risk_func a function which sums a vector across age and risk groups (created automatically in the code)
#' @param travel_matrix a square matrix with side length n_countries. Travel_matrix[x,y] is the proportion of time an individual in location x spends in location y
#' @param vax_allocation_params list of parameters for vaccine allocation
#' @param S numeric vector of length n. number of unvaccinated susceptibles.
#' @param n is number of locations * number of age groups * number of risk groups
#' @param E numeric vector of length n. number of unvaccinated exposed
#' @param I numeric vector of length n. number of unvaccinated infectious
#' @param R numeric vector of length n. number of unvaccinated recovered
#' @param vax_pool numeric vector of length 1. number of vaccines available 
#' @return n_vax_allocated numeric vector of length n_countries.  number of vaccines allocated to each country
#' @export
vaccinate_by_current_seasonal_alloc <- function(sum_age_risk_func,
                                          travel_matrix,
                                          vax_allocation_params,
                                          S, E, I, R, vax_pool) {
  # allocate proportional to seasonal coverage
  n_vax_allocated <- vax_allocation_params$coverage * vax_pool
  return(n_vax_allocated)
}

#' linear function of vaccine production
#' 
#' no vaccine produced until time vax_production_params[["detection_delay"]] + vax_production_params[["production_delay"]], then constant production rate until max number of doses ever made reached, then no production
#' @param vax_production_params list of parameters for vaccine production
#' @param t scalar time
#' @return needs to return a scalar: the number of vaccines ever produced up to time t
#' @export
produce_vax_linear_with_delay <- function(vax_production_params, t) {
  t_since_production <- t - (vax_production_params[["detection_delay"]] + 
                               vax_production_params[["production_delay"]])
  if(t_since_production < 0) {
    0
  } else {
    min(vax_production_params[["max_vax"]],
        t_since_production * vax_production_params[["production_rate"]])
  }
  }

######################################################################
# Functions for manipulating simulation output
######################################################################

#' @export
# this is tend... I presume the final time?
time_end <- function(results){
  ncol(results$S)
}

#' @export
# vector of deaths... What does each entry corresponds to? I presume split by
#  coutry, risk group, age group?
deaths <- function(results, popns){
  tend <- time_end(results)
  popns - results$S[,tend] - results$SV[,tend] - results$E[,tend] - results$EV[,tend] -
    results$I[,tend] - results$IV[,tend] - results$R[,tend] - results$RV[,tend]
}

#' @export
# total number of deaths
worldwide_deaths <- function(results, popns){
  sum(deaths(results, popns))
}

#' @export
# global attack rate
global_attack <- function(results, popns){
  pop_total <- sum(popns)
  tend <- time_end(results)
  sum(results$R[ ,tend] + results$RV[ ,tend] + deaths(results, popns)) / pop_total
}

#' @export
# Create a data frame from of the simulation results - worldwide deaths, global 
#  attack rate
deaths_GAR_df <- function(results, popns){
    data.frame("worldwide_deaths" = vapply(results, worldwide_deaths, double(1), popns), 
             "global_attack" = vapply(results, global_attack, double(1), popns))
  }


# attack rate by country
country_attack <- function(results, popns, labels){
  pop_total <- sum(popns)
  tend <- time_end(results)
  sum_age_risk <- sum_age_risk_closure(labels)
  final_size_by_group <- results$R[ ,tend] + results$RV[ ,tend] + deaths(results, popns)
  country_attack_rate <- sum_age_risk(final_size_by_group) / sum_age_risk(popns)
}

#' @export
# extract the number of deaths from each list element in the results list.
many_dead <- function(results, popns){
  dead <- lapply(results, deaths, popns = popns)
  dead <- do.call("cbind", dead)
  dead <- cbind(labels, dead)
}


