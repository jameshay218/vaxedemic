#' return the full results of main_simulation
#' 
#' @param res output of main_simulation
#' @return output of main_simulation
#' @export
return_all_res <- function(res){
  return(res)
}


#' @export
calc_region_time_series <- function(res, X, labels){
    time_series <- data.table::data.table(res$I)
    I <- combine_incidence_region(I, labels)
    tmp <- unique(I[,c("Location","variable","sumI","sumN")])
    return(tmp)
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
  return(list(peakTimes = peakTimes$peakTimes, 
              country_attack_rate = country_attack_rate$country_attack_rate))
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
