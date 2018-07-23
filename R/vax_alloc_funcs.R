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
#' n is number of locations * number of age groups * number of risk groups
#' @param E numeric vector of length n. number of unvaccinated exposed
#' @param I numeric vector of length n. number of unvaccinated infectious
#' @param R numeric vector of length n. number of unvaccinated recovered
#' @param SV numeric vector of length n. number of vaccinated susceptibles.
#' @param EV numeric vector of length n. number of vaccinated exposed
#' @param IV numeric vector of length n. number of vaccinated infectious
#' @param RV numeric vector of length n. number of vaccinated recovered
#' @param incidence numeric vector of length n. incidence
#' @param vax_pool numeric vector of length 1. number of vaccines available 
#' @return n_vax_allocated numeric vector of length n_countries.  number of vaccines allocated to each country
#' @export
vaccinate_by_incidence <- function(sum_age_risk_func, 
                                   travel_matrix,
                                   vax_allocation_params,
                                   S, E, I, R, 
                                   SV, EV, IV, RV, incidence, vax_pool) {
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
#' n is number of locations * number of age groups * number of risk groups
#' @param E numeric vector of length n. number of unvaccinated exposed
#' @param I numeric vector of length n. number of unvaccinated infectious
#' @param R numeric vector of length n. number of unvaccinated recovered
#' @param SV numeric vector of length n. number of vaccinated susceptibles.
#' @param EV numeric vector of length n. number of vaccinated exposed
#' @param IV numeric vector of length n. number of vaccinated infectious
#' @param RV numeric vector of length n. number of vaccinated recovered
#' @param incidence numeric vector of length n. incidence
#' @param vax_pool numeric vector of length 1. number of vaccines available 
#' @return n_vax_allocated numeric vector of length n_countries.  number of vaccines allocated to each country
#' @export
vaccinate_by_current_seasonal_alloc <- function(sum_age_risk_func,
                                                travel_matrix,
                                                vax_allocation_params,
                                                S, E, I, R, 
                                                SV, EV, IV, RV, incidence, vax_pool) {
  # allocate proportional to seasonal coverage
  n_vax_allocated <- vax_allocation_params$coverage * vax_pool
  return(n_vax_allocated)
}