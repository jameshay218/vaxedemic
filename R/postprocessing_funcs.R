#' postprocess simulations to plot the peak time and save outputs
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @return NULL
#' @export
postprocessing_peak_times_and_attack_rates <- function(res_list, runName, other_info) {
  # postprocessing: plot peak times
  
  # res_list$res is a list of length n_runs.  Each element is a list containing
  # two elements: peakTimes and country_attack_rate.
  # this lapply pulls out the peakTimes element for each run.
  res <- lapply(res_list$res, function(x) x$peakTimes)
  labels <- res_list$processed_inputs[["labels"]]
  regionDat <- other_info[["regionDat"]]
  latitudeDat <- other_info[["latitudeDat"]]
  
  median_ci <- calc_median_ci_by_country(res, labels, 
                                   regionDat, 
                                   latitudeDat)
  p <- plot_peak_times(median_ci)
  
  # save the plot
  filename <- paste0("outputs/peakTimes",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  # median and cis
  write.table(median_ci, paste0(filename, "_median_ci.csv"),sep=",",row.names=FALSE)
  
  # peak times for all runs
  res <- neaten_raw_output_by_country(res, labels, 
                                            regionDat, 
                                            latitudeDat)
  write.table(res, paste0(filename, "_all_runs.csv"),sep=",",row.names=FALSE)
  
  # postprocessing: plot attack rates
  res <- lapply(res_list$res, function(x) x$country_attack_rate)

  median_ci <- calc_median_ci_by_country(res, labels, 
                                   regionDat, 
                                   latitudeDat)
  p <- plot_country_attack_rates(median_ci)
  
  # save the plot
  filename <- paste0("outputs/country_attack_rates",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  write.table(median_ci, paste0(filename, "_median_ci.csv"),sep=",",row.names=FALSE)
  
  res <- neaten_raw_output_by_country(res, labels, 
                                      regionDat, 
                                      latitudeDat)
  write.table(res, paste0(filename, "_all_runs.csv"),sep=",",row.names=FALSE)
  return(NULL)
}

#' postprocess simulations to plot the peak time and save outputs
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @return NULL
#' @export
postprocessing_peak_times <- function(res_list, runName, other_info) {
  # postprocessing: plot peak times
  res <- res_list$res
  labels <- res_list$processed_inputs[["labels"]]
  regionDat <- other_info[["regionDat"]]
  latitudeDat <- other_info[["latitudeDat"]]
  
  median_ci <- calc_median_ci_by_country(res, labels, 
                                   regionDat, 
                                   latitudeDat)
  p <- plot_peak_times(median_ci)
  
  # save the plot
  filename <- paste0("outputs/",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  write.table(median_ci, paste0(filename, "_median_ci.csv"),sep=",",row.names=FALSE)
  
  res <- neaten_raw_output_by_country(res, labels, 
                                      regionDat, 
                                      latitudeDat)
  write.table(res, paste0(filename, "_all_runs.csv"),sep=",",row.names=FALSE)
  return(NULL)
}

#' postprocess simulations to plot the country attack rate and save outputs
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @param other_info list which provides any other information needed, such as to calculate the summaries
# or post-process results.
#' @return NULL
#' @export
postprocessing_country_attack <- function(res_list, runName, other_info) {
  # postprocessing: plot attack rates
  res <- res_list$res
  labels <- res_list$processed_inputs[["labels"]]
  regionDat <- other_info[["regionDat"]]
  latitudeDat <- other_info[["latitudeDat"]]
  
  median_ci <- calc_median_ci_by_country(res, labels, 
                                   regionDat, 
                                   latitudeDat)
  p <- plot_country_attack_rates(median_ci)
  
  # save the plot
  filename <- paste0("outputs/",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  write.table(median_ci, paste0(filename, "_median_ci.csv"),sep=",",row.names=FALSE)
  
  res <- neaten_raw_output_by_country(res, labels, 
                                      regionDat, 
                                      latitudeDat)
  write.table(res, paste0(filename, "_all_runs.csv"),sep=",",row.names=FALSE)
  return(NULL)
}

#' just save summary statistics from simulations as an .rds file
#' 
#' @param res_list summary statistics for each run and processed inputs produced by run_fixed_params
#' @param runName character string to make filename out of
#' @return NULL
#' @export
postprocessing_simple_save <- function(res_list, runName) {
  # postprocessing: plot peak times
  res <- res_list$res
  
  # save the plot
  filename <- paste0("outputs/calibration_",runName, "_output.rds")
  saveRDS(res, filename)
  return(NULL)
}