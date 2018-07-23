#' return the full results of main_simulation
#' 
#' @param res output of main_simulation
#' @return output of main_simulation
#' @export
return_all_res <- function(res){
  return(res)
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
  V_class <- c("SV", "EV", "IV", "RV")
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
  I <- data.table::data.table(res$I)
  I <- combine_incidence(I, labels)
  tmp <- unique(I[,c("Location","variable","sumI","sumN")])
  peakTimes <- plyr::ddply(tmp,~Location, function(x) x$variable[which.max(x$sumI)])[,2]
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
  pop_total <- sum(X)
  tend <- time_end(res)
  sum_age_risk <- sum_age_risk_closure(labels)
  final_size_by_group <- res$R[ ,tend] + res$RV[ ,tend] + deaths(res, X)
  country_attack_rate <- sum_age_risk(final_size_by_group) / sum_age_risk(X)
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
