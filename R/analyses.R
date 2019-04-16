######################################################################
# Functions for manipulating simulation output
######################################################################

#' @export
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

shrink_data <- function(res){
  compartments <- names(res)
  final <- NULL
  for(i in 1:length(res)){
    tmp <- data.table(res[[i]])
    tmp <- melt(tmp)
    tmp$comparment <- compartments[i]
    final[[i]] <- tmp
  }
  final <- rbindlist(final)
  return(final)    
}

combine_incidence <- function(I, labels){
    I <- cbind(labels, I)
    I <- melt(I, id.vars=colnames(labels))
    I$variable <- as.numeric(as.character(I$variable))
    I <- data.table(I)
    I[,sumI:=sum(value),key=c("Location","variable")]
    I[,sumN:=sum(X),key=c("Location","variable")]
    return(I)
}

combine_incidence_region <- function(I, labels){
    I <- cbind(labels, I)
    I <- melt(I, id.vars=colnames(labels))
    I$variable <- as.numeric(as.character(I$variable))
    I <- data.table(I)
    I[,sumI:=sum(value),key=c("region","variable")]
    I[,sumN:=sum(X),key=c("region","variable")]
    return(I)
}

#' @export
# Create a data frame from of the simulation results - worldwide deaths, global 
#  attack rate
deaths_GAR_df <- function(results, popns){
  data.frame("worldwide_deaths" = vapply(results, worldwide_deaths, double(1), popns), 
             "global_attack" = vapply(results, global_attack, double(1), popns))
}

#' @export
# extract the number of deaths from each list element in the results list.
many_dead <- function(results, popns){
  dead <- lapply(results, deaths, popns = popns)
  dead <- do.call("cbind", dead)
  dead <- cbind(labels, dead)
}

#' calculate the mean, median and 2.5%, 97.5% percentiles of a summary statistic by country
#' 
#' @param res summary statistics outputted by run_simulation.  Should be a list
#' of length n_runs, where each element is a numeric vector of length n_countries.
#' @param labels data frame outputted by setup_inputs containing country names
#' @param regionDat data frame containing a column for Location and a column for
#' the region that location is in. read from data/regions_clean.csv
#' @param latitudeDat data frame containing a column for Location and a column for
#' the latitude of the location. read from data/latitudes_intersect.csv
#' @return a data frame with the mean, median, 2.5% and 97.5% percentiles for
#' each country; the region and hemisphere each country belongs to; and the latitude
#' of each country
calc_median_ci_by_country <- function(res, labels, regionDat, latitudeDat) {
  res = do.call("cbind",res)
  summary_stats <- t(apply(res, 1, function(x) c(mean(x),quantile(x, c(0.025,0.5,0.975)))))
  colnames(summary_stats) <- c("mean","lower95","median","upper95")
  summary_stats <- data.frame(Location=unique(labels$Location), summary_stats)
  dat <- merge(summary_stats,regionDat[,c("Location","region")])
  dat <- merge(dat,latitudeDat)
  tmp <- unique(dat[,c("Location","latitude")])
  my_latitudes <- order(tmp$latitude)
  dat$Location <- factor(dat$Location, levels=tmp$Location[my_latitudes])
  dat$Hemisphere <- ifelse(dat$latitude < 0,"Southern","Northern")
  return(dat)
}

#' format the summary statistic by country over all n_runs
#' 
#' @param res summary statistics outputted by run_simulation.  Should be a list
#' of length n_runs, where each element is a numeric vector of length n_countries.
#' @param labels data frame outputted by setup_inputs containing country names
#' @param regionDat data frame containing a column for Location and a column for
#' the region that location is in. read from data/regions_clean.csv
#' @param latitudeDat data frame containing a column for Location and a column for
#' the latitude of the location. read from data/latitudes_intersect.csv
#' @param melted if TRUE, returns the data frame melted by the outcome variable of interest.
#' otherwise, enumerates outcome variable across columns
#' @return a data frame with the summary statistic for each run for
#' each country; the region and hemisphere each country belongs to; and the latitude
#' of each country
neaten_raw_output_by_country <- function(res, labels, regionDat, latitudeDat, melted=FALSE) {
  res = do.call("cbind",res)
  colnames(res) <- seq_len(ncol(res))
  summary_stats <- data.frame(Location=unique(labels$Location), res)
  dat <- merge(summary_stats,regionDat[,c("Location","region")])
  dat <- merge(dat,latitudeDat)
  tmp <- unique(dat[,c("Location","latitude")])
  my_latitudes <- order(tmp$latitude)
  dat$Location <- factor(dat$Location, levels=tmp$Location[my_latitudes])
  dat$Hemisphere <- ifelse(dat$latitude < 0,"Southern","Northern")
  if(melted){
      dat <- reshape2::melt(dat, id.vars=setdiff(colnames(dat),paste0("X",colnames(res))))
      dat$Location <- factor(dat$Location, levels=tmp$Location[my_latitudes])
  }
  
  return(dat)
}

calc_median_global_deaths <- function(filename) {
  median(calc_global_deaths(filename))
}

calc_global_deaths <- function(filename) {
  deaths <- readRDS(filename)
  vnapply(deaths, function(x) sum(x[, ncol(x)]))
}

