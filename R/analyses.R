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

######################################################################
                                        # Functions for manipulating simulation output
######################################################################

                                        # this is tend... I presume the final time?
#' @export
time_end <- function(results){
    ncol(results$S)
}


                                        # vector of deaths... What does each entry corresponds to? I presume split by
                                        #  coutry, risk group, age group?
#' @export
deaths <- function(results){
    tend <- time_end(results)
    popns - results$S[,tend] - results$SV[,tend] - results$E[,tend] - results$EV[,tend] -
        results$I[,tend] - results$IV[,tend] - results$R[,tend] - results$RV[,tend]
}

                                        # total number of deaths
#' @export
worldwide_deaths <- function(results){
    sum(deaths(results))
}

                                        # global attack rate
#' @export
global_attack <- function(results){
    pop_total <- sum(popns)
    tend <- time_end(results)
    sum(results$R[,tend] + results$RV[,tend] + deaths(results))/pop_total
}

#' @export
find_peak_times_list <- function(res){
    return(lapply(res, find_peak_times))
    
}

#' @export
find_peak_times <- function(res){
    peakTimes <- (apply(res$I, 1, function(x) as.numeric(colnames(res$I)[which.max(x)])))
    return(peakTimes)
}

find_peak_summaries <- function(res, labels){
    peakTimes <- find_peak_times_list(res)
    peakTimes <- do.call("cbind",peakTimes)
    summaryPeaks <- t(apply(peakTimes, 1, function(x) c(mean(x),quantile(x, c(0.025,0.5,0.975)))))
    colnames(summaryPeaks) <- c("mean","lower95","median","upper95")
    summaryPeaks <- cbind(labels, summaryPeaks)
    return(summaryPeaks)
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

calculate_summaries <- function(res, labels, requested_stats,...){
    if(requested_stats == "all_res"){
        return(res)
    }    
    I <- data.table(res$I)
    I <- combine_incidence(I, labels)
    tmp <- unique(I[,c("Location","variable","sumI","sumN")])
    peakTimes <- ddply(tmp,~Location, function(x) x$variable[which.max(x$sumI)])[,2]
    return(peakTimes)
    #tmp <- merge(tmp, regionDat[,c("Location","region")], by="Location")
    tmp <- merge(tmp, latitudeDat[,c("Location","latitude")],by="Location")
    tmp$loc <- ifelse(tmp$latitude > 0, "North", "South")
    tmp[abs(tmp$latitude) < 23.5,"loc"] <- "Tropics"
    tmp[,I:=sum(sumI),key=c("loc","variable")]
    tmp[,N:=sum(sumN),key=c("loc","variable")]
    tmp <- unique(tmp[,c("variable","I","N","loc")])    
    return(list("I"=tmp,"peakTimes"=peakTimes))

}

                                        # Create a data frame from of the simulation results - worldwide deaths, global 
                                        #  attack rate
deaths_GAR_df <- function(results){
    data.frame("worldwide_deaths" = vapply(results, worldwide_deaths, double(1)), 
               "global_attack" = vapply(results, global_attack, double(1)))
}


