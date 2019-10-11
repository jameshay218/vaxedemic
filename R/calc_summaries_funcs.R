#' return the full results of main_simulation
#' 
#' @param res output of main_simulation
#' @return output of main_simulation
#' @export
return_all_res <- function(res){
  return(res)
}

#' calculate the incidence, peak times, attack rates and deaths by country
#'
#' @param res output of main_simulation
#' @param X the vector of population sizes in each location/age/risk group
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @return a list with four entries.
#' incidence: the incidence over time for each country
#' vaccinated: the number of vaccinated individuals over time for each country
#' peakTimes: numeric vector of length n_countries: the peak time for each country
#' country_attack_rate: numeric vector of length n_countries: the attack rate for each country
#' @export
calc_incidence_peak_times_attack_rates_deaths <- function(res, X, labels){
  incidence <- calc_incidence_time_series(res, X, labels)
  peakTimes <- calc_peak_times(res, labels)
  country_attack_rate <- calc_country_attack(res, X, labels)
  deaths <- calc_cumulative_deaths_time_series(res, X, labels)
  return(c(incidence, peakTimes, country_attack_rate, deaths))
}

#' calculate the incidence, number of vaccinated individuals, peak times
#' and attack rates by country
#'
#' @param res output of main_simulation
#' @param X the vector of population sizes in each location/age/risk group
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @return a list with four entries.
#' incidence: the incidence over time for each country
#' vaccinated: the number of vaccinated individuals over time for each country
#' peakTimes: numeric vector of length n_countries: the peak time for each country
#' country_attack_rate: numeric vector of length n_countries: the attack rate for each country
#' @export
calc_incidence_vaccinated_peak_times_attack_rates <- function(res, X, labels){
  incidence <- calc_incidence_time_series(res, X, labels)
  vaccinated <- calc_vaccinated_time_series(res, X, labels)
  peakTimes <- calc_peak_times(res, labels)
  country_attack_rate <- calc_country_attack(res, X, labels)
  return(c(incidence, vaccinated, peakTimes, country_attack_rate))
}

#' calculate the incidence, number of vaccinated individuals, peak times,
#' attack rates and deaths by country
#'
#' @param res output of main_simulation
#' @param X the vector of population sizes in each location/age/risk group
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @return a list with four entries.
#' incidence: the incidence over time for each country
#' vaccinated: the number of vaccinated individuals over time for each country
#' peakTimes: numeric vector of length n_countries: the peak time for each country
#' country_attack_rate: numeric vector of length n_countries: the attack rate for each country
#' @export
calc_incidence_vaccinated_peak_times_attack_rates_deaths <- function(res, X, labels){
  incidence <- calc_incidence_time_series(res, X, labels)
  vaccinated <- calc_vaccinated_time_series(res, X, labels)
  peakTimes <- calc_peak_times(res, labels)
  country_attack_rate <- calc_country_attack(res, X, labels)
  deaths <- calc_cumulative_deaths_time_series(res, X, labels)
  return(c(incidence, vaccinated, peakTimes, country_attack_rate, deaths))
}

#' calculate the peak time and attack rate by country
#'
#' @param res output of main_simulation
#' @param X the vector of population sizes in each location/age/risk group
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @return a list with two entries.
#' peakTimes: numeric vector of length n_countries: the peak time for each country
#' country_attack_rate: numeric vector of length n_countries: the attack rate for each country
#' @export
calc_peak_times_and_attack_rates <- function(res, X, labels){
  peakTimes <- calc_peak_times(res, labels)
  country_attack_rate <- calc_country_attack(res, X, labels)
  return(c(peakTimes, country_attack_rate))
}

#' calculate the incidence over time for each country
#'
#' @param res output of main_simulation
#' @param X the vector of population sizes in each location/age/risk group
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @return a list containing a matrix with n_countries rows. 
#' Each row contains the incidence over time for that country.
#' @export
calc_incidence_time_series <- function(res, X, labels){
  incidence <- calc_country_time_series(res, X, labels, "incidence")
  return(list(incidence = incidence))
}

#' calculate the number of vaccinated individuals over time for each country
#'
#' @param res output of main_simulation
#' @param X the vector of population sizes in each location/age/risk group
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @return a list containing a matrix with n_countries rows. 
#' Each row contains the incidence over time for that country.
#' @export
calc_vaccinated_time_series <- function(res, X, labels){
  V_class <- c("SV1", "EV1", "IV1", "RV1", "SV2", "EV2", "IV2", "RV2",
               "SP", "EP", "IP", "RP")
  vaccinated <- calc_country_time_series(res, X, labels, V_class)
  return(list(vaccinated = vaccinated))
}

#' calculate the peak time by country
#'
#' @param res output of main_simulation
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @return a list containing a numeric vector of length n_countries: 
#' the peak time for each country
#' @export
calc_peak_times <- function(res, labels){
  sum_age_risk <- sum_age_risk_closure(labels)
  country_incidence <- apply(res$incidence, 2, sum_age_risk)
  daily_incidence <- thin_time_series(country_incidence, 
                                                   thin_integer = TRUE,
                                                   thin_by_sum = TRUE)
  weekly_incidence <- thin_time_series(daily_incidence, 
                                      thin_every = 7,
                                      thin_by_sum = TRUE)
  peakTimes_idx <- apply(weekly_incidence, 1, which.max)
  time_vec <- as.numeric(colnames(weekly_incidence))
  peakTimes <- time_vec[peakTimes_idx]
  return(list(peakTimes = peakTimes))
}

#' calculate the attack rate by country
#'
#' @param res output of main_simulation
#' @param X the vector of population sizes in each location/age/risk group
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @return a list ocntaining a numeric vector of length n_countries: 
#' the attack rate for each country
#' @export
calc_country_attack <- function(res, X, labels){
  incidence <- calc_country_time_series(res, X, labels, "incidence")
  final_size <- rowSums(incidence)
  pop_total <- sum(X)
  tend <- time_end(res)
  sum_age_risk <- sum_age_risk_closure(labels)
  country_attack_rate <- final_size / sum_age_risk(X)
  return(list(country_attack_rate = country_attack_rate))
}

#' calculate the number of individuals summed over compartments over time for each country
#'
#' @param res output of main_simulation
#' @param X the vector of population sizes in each location/age/risk group
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @param compartments vector of compartment names to sum over
#' @return a list containing a matrix with n_countries rows. 
#' Each row contains the number of individuals summed over compartments over time for that country.
calc_country_time_series <- function(res, X, labels, compartments){
  pop_total <- sum(X)
  tend <- time_end(res)
  sum_age_risk <- sum_age_risk_closure(labels)
  
  if(length(compartments) > 1) {
    res <- sum_list(res[compartments])
  } else {
    res <- res[[compartments]]
  }
  
  country_time_series <- apply(res, 2, sum_age_risk)
  rownames(country_time_series) <- levels(labels$Location)
  return(country_time_series)
}

#' calculate the cumulative number of deaths over time for each country
#'
#' @param res output of main_simulation
#' @param X the vector of population sizes in each location/age/risk group
#' @param labels data frame outputted by setup_inputs containing the location/
#' age/risk group for each row in res
#' @return a list containing a matrix with n_countries rows. 
#' Each row contains the incidence over time for that country.
#' @export
calc_cumulative_deaths_time_series <- function(res, X, labels){
  deaths <- calc_country_time_series(res, X, labels, "deaths")
  cumulative_deaths <- t(apply(deaths, 1, cumsum))
  return(list(dead = cumulative_deaths))
}
