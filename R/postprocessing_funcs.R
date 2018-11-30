#' postprocess simulations to plot the peak time and attack rates,
#' and save outputs for the incidence, vaccinated, peak time and attack rates
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @param output_prefix character vector of length 1.  Prefix for output filenames
#' @return NULL
#' @export
postprocessing_incidence_vaccinated_peak_times_attack_rates <- function(res_list, runName, other_info, output_prefix) {
  postprocessing_save_country_incidence(res_list, runName, other_info, output_prefix)
  postprocessing_save_country_vaccinated(res_list, runName, other_info, output_prefix)
  postprocessing_peak_times(res_list, runName, other_info, output_prefix)
  postprocessing_country_attack(res_list, runName, other_info, output_prefix)
  return(NULL)
}

#' postprocess simulations to plot the peak time and attack rates, and save outputs
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @param output_prefix character vector of length 1.  Prefix for output filenames
#' @return NULL
#' @export
postprocessing_peak_times_and_attack_rates <- function(res_list, runName, other_info, output_prefix) {
  postprocessing_peak_times(res_list, runName, other_info, output_prefix)
  postprocessing_country_attack(res_list, runName, other_info, output_prefix)
  return(NULL)
}

#' postprocess simulations to save the daily incidence for each country
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @param output_prefix character vector of length 1.  Prefix for output filenames
#' @return NULL
#' @export
postprocessing_save_country_incidence <- function(res_list, runName, other_info, output_prefix){
  
  # thin country times series to one observation per day before saving
  thin_incidence <- function(incidence) {
    thin_time_series(incidence, thin_integer = TRUE, thin_by_sum = TRUE)
  }
  
  res_list$res <- lapply(res_list$res, function(x) thin_incidence(x$incidence))
  # save country time series as rds
  saveRDS(res_list$res, paste0(output_prefix, "_",runName,"_incidence.rds"))
  return(NULL)
}

#' postprocess simulations to save the daily number of vaccinated individuals for each country
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @param output_prefix character vector of length 1.  Prefix for output filenames
#' @return NULL
#' @export
postprocessing_save_country_vaccinated <- function(res_list, runName, other_info, output_prefix){
  
  # thin country times series to one observation per day before saving
  thin_vaccinated <- function(vaccinated) {
    thin_time_series(vaccinated, thin_integer = TRUE, thin_by_sum = FALSE)
  }
  
  res_list$res <- lapply(res_list$res, function(x) thin_vaccinated(x$vaccinated))
  # save country time series as rds
  saveRDS(res_list$res, paste0(output_prefix, "_",runName,"_vaccinated.rds"))
  return(NULL)
}

#' postprocess simulations to plot the peak time and save outputs
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @param output_prefix character vector of length 1.  Prefix for output filenames
#' @return NULL
#' @export
postprocessing_peak_times <- function(res_list, runName, other_info, output_prefix) {
  # postprocessing: plot peak times
  
  # res_list$res is a list of length n_runs.  Each element is a list containing
  # two elements: peakTimes and country_attack_rate.
  # this lapply pulls out the peakTimes element for each run.
  res_PT <- lapply(res_list$res, function(x) x$peakTimes)
  labels <- res_list$processed_inputs[["labels"]]
  regionDat <- other_info[["regionDat"]]
  latitudeDat <- other_info[["latitudeDat"]]
  
  median_ci <- calc_median_ci_by_country(res_PT, labels, 
                                         regionDat, 
                                         latitudeDat)
  p <- plot_peak_times(median_ci)
  
  # save the plot
  filename <- paste0(output_prefix, "_peakTimes_",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  # median and cis
  write.table(median_ci, paste0(filename, "_median_ci.csv"),sep=",",row.names=FALSE)
  
  # peak times for all runs
  res_PT_clean<- neaten_raw_output_by_country(res_PT, labels, 
                                              regionDat, 
                                              latitudeDat)
  write.table(res_PT_clean, paste0(filename, "_all_runs.csv"),sep=",",row.names=FALSE)
  
  ## melted, for density plot
  ## To plot density strips of attack rates
  res_PT_melted<- neaten_raw_output_by_country(res_PT, labels, 
                                               regionDat, 
                                               latitudeDat,
                                               melt=TRUE)
  ## If we want to save the melted raw output - commented out as this is a large file
  #write.table(res_PT_melted, paste0(filename, "_all_runs_melted.csv"),sep=",",row.names=FALSE)
  p_PT_dens <- density_country_peak_times(res_PT_melted)
  # save the plot
  filename <- paste0(output_prefix, "_peakTimes_density_",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p_PT_dens)
  dev.off()
  return(NULL)
}

#' postprocess simulations to plot the country attack rate and save outputs
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @param output_prefix character vector of length 1.  Prefix for output filenames
#' @return NULL
#' @export
postprocessing_country_attack <- function(res_list, runName, other_info, output_prefix) {
  # postprocessing: plot attack rates
  res_AR <- lapply(res_list$res, function(x) x$country_attack_rate)
  labels <- res_list$processed_inputs[["labels"]]
  regionDat <- other_info[["regionDat"]]
  latitudeDat <- other_info[["latitudeDat"]]
  
  median_ci <- calc_median_ci_by_country(res_AR, labels, 
                                         regionDat, 
                                         latitudeDat)
  p <- plot_country_attack_rates(median_ci)
  
  # save the plot
  filename <- paste0(output_prefix, "_country_attack_rates_",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  write.table(median_ci, paste0(filename, "_median_ci.csv"),sep=",",row.names=FALSE)
  
  res_AR_clean<- neaten_raw_output_by_country(res_AR, labels, 
                                              regionDat, 
                                              latitudeDat)    
  write.table(res_AR_clean, paste0(filename, "_all_runs.csv"),sep=",",row.names=FALSE)
  
  ## To plot density strips of peak times
  res_AR_melted<- neaten_raw_output_by_country(res_AR, labels, 
                                               regionDat, 
                                               latitudeDat,
                                               melt=TRUE)
  ## If we want to save the melted raw output - commented out as this is a large file
  #write.table(res_AR_melted, paste0(filename, "_all_runs_melted.csv"),sep=",",row.names=FALSE)
  p_AR_dens <- density_country_attack_rates(res_AR_melted)
  # save the plot
  filename <- paste0(output_prefix, "_country_attack_rates_density_",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p_AR_dens)
  dev.off() 
  return(NULL)
}

#' just save summary statistics from simulations as an .rds file
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param output_prefix character vector of length 1.  Prefix for output filenames
#' @return NULL
#' @export
postprocessing_simple_save <- function(res_list, runName) {

  # postprocessing: plot peak times
  res <- res_list$res
  
                                        # save the plot
    filename <- paste0(output_prefix, "_",runName, "_output.rds")
    filename <- paste0(output_prefix, "_",runName, "_output.RData")
    save(res_list, file=filename)
                                        #  saveRDS(res, filename)
  return(NULL)
}
