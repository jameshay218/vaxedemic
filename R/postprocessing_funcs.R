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
  
  dat <- calc_median_ci_by_country(res, labels, 
                                   regionDat, 
                                   latitudeDat)
  p <- plot_peak_times(dat)
  
  # save the plot
  filename <- paste0("outputs/calibration_",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  write.table(dat, paste0(filename, "_data.csv"),sep=",",row.names=FALSE)
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
  # postprocessing: plot peak times
  res <- res_list$res
  labels <- res_list$processed_inputs[["labels"]]
  regionDat <- other_info[["regionDat"]]
  latitudeDat <- other_info[["latitudeDat"]]
  
  dat <- calc_median_ci_by_country(res, labels, 
                                   regionDat, 
                                   latitudeDat)
  p <- plot_country_attack_rates(dat)
  
  # save the plot
  filename <- paste0("outputs/calibration_",runName)
  
  png(paste0(filename, "_plot.png"),width=800,height=1200)
  plot(p)
  dev.off()
  
  # save the output
  write.table(dat, paste0(filename, "_data.csv"),sep=",",row.names=FALSE)
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