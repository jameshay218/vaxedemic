#' Simple plot
#'
#' Simple plot of simulation outputs
#' @import data.table
#' @export
model_plot_simple <- function(res, labels, n_countries){
  times <- as.numeric(colnames(res$I))
  I <- cbind(labels[,c("Location","Age","RiskGroup")], res$I + res$IV)
  I <- melt(I, id.vars=c("Location","Age","RiskGroup"))
  I$Age <- as.factor(I$Age)
  I$RiskGroup <- as.factor(I$RiskGroup)
  I$variable <- as.numeric(I$variable)
  I$variable <- times[I$variable]
  I <- data.table(I)
  I_aggregated <- I[,x:=sum(value),by=list(variable,Location,Age)]
  N <- data.table(aggregate(data=labels, X~Location + Age,FUN=sum))
  N <- N[N$Location %in% unique(I_aggregated$Location),]
  I_aggregated <- merge(I_aggregated,N,by=c("Location","Age"))
  use_countries <- sample(unique(I_aggregated$Location), n_countries)
  #    ggplot(I_aggregated[I_aggregated$Location %in% unique(I_aggregated$Location)[n_countries],],aes(x=variable,y=x/X,col=Age)) +
  ggplot(I_aggregated[I_aggregated$Location %in% use_countries,],aes(x=variable,y=x/X,col=Age)) +
    geom_line() +
    facet_wrap(~Location) +
    theme_bw()
}


#' Model visualisations - MAIN FUNCTION
#'
#' TO DO: 1) need to make sure that the denominator is correct and that we are grouping correctly; 
#' 2) need to generalise grouping indices
#' 3) Make each country occupy same space in y-axis in heatmap plot
#' This is very much a work in progress. Takes the output of \code{\link{run_simulation}} and returns 
#' a set of incidence curve and heatmaps for per capita incidence by region. Requires data table and ggplot2.
#' @param the list of simulation outputs as in \code{\link{run_simulation}}
#' @param n_ages the number of age groups
#' @param n_riskgroups the number of risk groups
#' @param countryFile the full file path to the clean country name file
#' @param latitudeFile the full file path to the clean country latitude file
#' @param regionFile the full file path to the clean country by region file
#' @param vitalFile the full file path to the clean country demography (ie. needs "N") file
#' @import data.table
#' @export
model_plots <- function(res, n_ages, n_riskgroups, 
                        countryFile="data/countries_intersect.csv", 
                        latitudeFile="data/latitudes_intersect.csv", 
                        regionFile="data/regions_clean.csv", 
                        vitalFile="data/demographic_data_intersect.csv"){
  ## This was from development
  #data <- readRDS("~/Documents/vaxedemic_tmp/results.rds")
  #list2env(data, globalenv())
  
  countries <- read.csv(countryFile,stringsAsFactors = FALSE,header=FALSE)
  n_countries <- nrow(countries)
  
  latitudeDat <- read.csv(latitudeFile,stringsAsFactors=FALSE)
  regionDat <- read.csv(regionFile,stringsAsFactors=FALSE)
  
  ## Vital dat should have the first column name as "Location" anyway, so
  ## might need to change this
  vitalDat <- read.csv(vitalFile,stringsAsFactors=FALSE)
  colnames(vitalDat)[colnames(vitalDat) == "countryID"] <- "Location"
  
  plot_labels <- expand.grid("Location"=unique(countries[,1]),"Age"=1:n_ages,"RiskGroup"=1:n_riskgroups)
  
  
  compartment_names <- c("S","E","I","R","SV","EV","IV","RV")
  all_popns <- list(res$S,res$E,res$I,res$R,res$SV,res$EV,res$IV,res$RV)
  all_popns <- lapply(all_popns, data.table)
  
  return(list(all_popns,plot_labels,compartment_names,regionDat,latitudeDat,vitalDat))
  
  melted_pops <- lapply(all_popns, function(x){
    x <- cbind(x, plot_labels)
    x <- melt(x, value.name="y",id=c("Location","Age","RiskGroup"),as.is=TRUE)
    x$variable <- as.numeric(as.character(x$variable))
    x$Age <- as.factor(x$Age)
    x$RiskGroup <- as.factor(x$RiskGroup)
    x
  })
  
  for(i in 1:length(melted_pops)) melted_pops[[i]][,compartment:=compartment_names[i]]
  popns <- rbindlist(melted_pops)
  popns <- data.table(popns)
  regionDat <- data.table(regionDat)
  latitudeDat <- data.table(latitudeDat)
  vitalDat <- data.table(vitalDat)
  
  popns <- merge(popns, regionDat[,c("Location","region")], by="Location")
  popns <- merge(popns, latitudeDat[,c("Location","latitude")],by="Location")
  popns <- merge(popns, vitalDat[,c("Location","N")],by="Location")
  popns <- popns[popns$Location %in% c("China","United_Kingdom","Australia"),]
  #ggplot(popns) + stat_summary(aes(x=variable,y=y,col=compartment),fun.y= sum, geom="line") + facet_wrap(~region, scales = "free_y")
  
  ## Group by IV and I, as well as by risk group and country
  region_popns <- as.data.frame(popns[popns$compartment %in% c("IV","I"),
                                      j=list(y=sum(y,na.rm=TRUE)),by=c("Age","region","variable")])
  tmp <- unique(popns[,c("region","N")])
  region_N <- as.data.frame(tmp[,j=list(regionN=sum(as.numeric(N))),by=region])
  region_popns <- merge(region_popns, region_N,by="region")
  
  p1 <- ggplot(region_popns) + 
    geom_line(aes(x=variable,y=y/regionN,col=Age)) + 
    facet_wrap(~region,nrow=1) + 
    theme_bw()
  
  
  by_latitude <- as.data.frame(popns[popns$compartment %in% c("IV","I"),
                                     j=list(y=sum(y,na.rm=TRUE)),by=c("Location","Age","variable","latitude","region","N")])
  by_latitude$propn <- by_latitude$y/by_latitude$N
  tmp <- unique(by_latitude[,c("Location","latitude")])
  my_latitudes <- order(tmp$latitude)
  by_latitude$Location <- factor(by_latitude$Location, levels=tmp$Location[my_latitudes])
  
  p2 <- ggplot(by_latitude,aes(x=variable,y=Location,fill=logP),stat="identity") + 
    geom_tile() + 
    scale_fill_gradient(low="blue",high="red") +
    facet_wrap(~region,scales="free_y", nrow=2)
  return(list(p1,p2))
}


#' @export
plot_age_group <- function(Imat, age_group, n, ymax=500){
  tmpImat <- Imat[seq(age_group,nrow(Imat),by=n),]
  tmpImat <- reshape2::melt(tmpImat)
  colnames(tmpImat) <- c("index","x","value")
  p <- ggplot(tmpImat) +
    geom_line(aes(x=x,y=value)) +
    facet_wrap(~index) +
    scale_y_continuous(limits=c(0,ymax)) +
    theme_bw()
  return(p)
}

#' @export
generate_data_labels <- function(tmax, tdiv, n_countries, n_ages, n_riskgroups){
  plot_labels <- expand.grid("Time"=seq(0,tmax,by=1/tdiv),"Location"=1:n_countries,"Age"=1:n_ages,"RiskGroup"=1:n_riskgroups)
}



format_data_for_plots <- function(data,
                                  labels,
                                  times,
                                  groupingCols=c("Location","Age","RiskGroup"),
                                  aggregateCols=c("Location","Age")
){
  dat <- cbind(labels[,groupingCols],data)
  meltedDat <- reshape2::melt(dat, id.vars=groupingCols)
  tmp <- apply(meltedDat[,groupingCols], 2, as.factor)
  I <- cbind(labels[,c("Location","Age","RiskGroup")], res$I)
  I <- melt(I, id.vars=c("Location","Age","RiskGroup"))
  
  if(!is.null(meltedDat$Age)){
    meltedDat$Age <- as.factor(meltedDat$Age)
  }
  if(!is.null(meltedDat$RiskGroup)){
    meltedDat$Age <- as.factor(meltedDat$RiskGroup)
  }
  if(!is.null(meltedDat$Location)){
    meltedDat$Age <- as.factor(meltedDat$Location)
  }
  
  if(!is.null(meltedDat$variable)){
    meltedDat$variable <- as.numeric(meltedDat$variable)
    meltedDat$variable <- times[meltedDat$variable]
  }
  
  dat_aggregated <- aggregate(dat[,"value"], dat[,c("variable",aggregateCols)], FUN=sum)
  N <- aggregate(data=labels)
  
  #I_aggregated <- aggregate(I[,"value"], I[,c("variable","Location","Age")], FUN=sum)
  
  N <- aggregate(data=labels, X~Location + Age,FUN=sum)
  I_aggregated <- merge(I_aggregated,N,id.vars=c("Location","Age"))
  
  p1 <- ggplot(I_aggregated,aes(x=variable,y=x/X,col=Age)) +
    geom_line() +
    facet_wrap(~Location) +
    theme_bw()
  
  p2 <- ggplot(I, aes(x=variable,y=value,col=RiskGroup)) + geom_line() + facet_grid(Age~Location) + theme_bw()
  
}

#' plot the median and 95% quantiles of the peak time in each country across simulations
#' 
#' @param dat a data frame outputted by calc_median_ci_by_country.
#' contains the mean, median, 2.5% and 97.5% percentiles for the peak time in
#' each country; the region and hemisphere each country belongs to; and the latitude
#' of each country
#' @return a ggplot object
#'@import ggplot2
plot_peak_times <- function(dat){
  ggplot(dat) +
    geom_errorbarh(aes(y=Location,x=median,xmax=upper95,xmin=lower95, col=Hemisphere)) +
    geom_point(aes(y=Location,x=median),size=0.5) +
    scale_x_continuous(limits=c(0,365)) +
    xlab("Peak time (days), median and 95% quantiles") +
    facet_grid(region~.,scales="free_y", space="free",switch="both") +
    theme(axis.text.y=element_text(size=6)) + theme_bw()
}

#' plot the median and 95% quantiles of the attack rate in each country across simulations
#' 
#' @param dat a data frame outputted by calc_median_ci_by_country.
#' contains the mean, median, 2.5% and 97.5% percentiles for the attack rate in
#' each country; the region and hemisphere each country belongs to; and the latitude
#' of each country
#' @return a ggplot object
#'@import ggplot2
plot_country_attack_rates <- function(dat) {
  ggplot(dat) +
    geom_errorbarh(aes(y=Location,x=median,xmax=upper95,xmin=lower95, col=Hemisphere)) +
    geom_point(aes(y=Location,x=median),size=0.5) +
    scale_x_continuous(limits=c(0,1)) +
    xlab("Attack rate, median and 95% quantiles") +
    facet_grid(region~.,scales="free_y", space="free",switch="both") +
    theme(axis.text.y=element_text(size=6)) + theme_bw()
}

#' plot the density strips of the attack rate in each country across simulations
#' 
#' @param dat a data frame outputted neaten_raw_output_by_country with
#' "melted=TRUE"
#' @return a ggplot object
#'@import ggplot2
density_country_attack_rates <- function(dat) {
  ggplot(dat, aes(x=value,y=factor(Location))) + 
    stat_density(aes(fill=..scaled..),geom="tile",position="identity")+ 
    scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
    scale_fill_gradient(low="white",high="black") +
    labs(fill="Scaled density") +
    ylab("Location") + 
    theme_bw() +
    theme(panel.grid=element_blank(),
          legend.position="bottom",
          legend.title.align = 0.5)+
    facet_grid(region~.,scales="free_y", space="free",switch="both")
  
}


#' plot the density strips of the peak times in each country across simulations
#' 
#' @param dat a data frame outputted neaten_raw_output_by_country with
#' "melted=TRUE"
#' @return a ggplot object
#'@import ggplot2
density_country_peak_times <- function(dat) {
  ggplot(dat, aes(x=value,y=factor(Location))) + 
    stat_density(aes(fill=..scaled..),geom="tile",position="identity")+
    scale_x_continuous(limits=c(0,365),expand=c(0,0)) +
    scale_fill_gradient(low="white",high="black") +
    labs(fill="Scaled density") +
    ylab("Location") + 
    theme_bw() +
    theme(panel.grid=element_blank(),
          legend.position="bottom",
          legend.title.align = 0.5)+
    facet_grid(region~.,scales="free_y", space="free",switch="both")
  
}

#' plot incidence time series for a country
#' 
#' plot 95% CI and ten runs
#' 
#' @param output list of n matrices, where n is the number of runs.
#' Each row in each matrix corresponds to a country and each column to a time.
#' the rownames of the matrices are the country names.
#' @param country_name name of country for which to plot time series
#' @param thin_every an integer.  If 0, ignore.  Otherwise, thin each matrix
#' every this many columns
#' @param thin_integer a logical. If TRUE, thin each matrix to integer values of 
#' times.  If FALSE, ignore
#' @return ggplot object with 95% CI and ten runs
#' @export
plot_country_incidence <- function(output, country_name,  
                                     thin_every = 7, 
                                     thin_integer = FALSE) {
  
  y_label <- make_y_label("incidence", thin_every, thin_integer)
  
  plot_country_time_series(output, 
                           country_name, 
                           y_label = y_label, 
                           thin_every = thin_every, 
                           thin_integer = thin_integer,
                           thin_by_sum = TRUE)
}

#' plot vaccinated individuals time series for a country
#' 
#' plot 95% CI and ten runs
#' 
#' @param output list of n matrices, where n is the number of runs.
#' Each row in each matrix corresponds to a country and each column to a time.
#' the rownames of the matrices are the country names.
#' @param country_name name of country for which to plot time series
#' @param thin_every an integer.  If 0, ignore.  Otherwise, thin each matrix
#' every this many columns
#' @param thin_integer a logical. If TRUE, thin each matrix to integer values of 
#' times.  If FALSE, ignore
#' @return ggplot object with 95% CI and ten runs
#' @export
plot_country_vaccinated <- function(output, country_name,  
                                   thin_every = 7, 
                                   thin_integer = FALSE) {

  y_label <- make_y_label("vaccinated", thin_every, thin_integer)
  
  plot_country_time_series(output, 
                           country_name, 
                           y_label = y_label, 
                           thin_every = thin_every, 
                           thin_integer = thin_integer,
                           thin_by_sum = FALSE)
}

#' plot time series for a country
#' 
#' plot 95% CI and ten runs
#' 
#' @param output list of n matrices, where n is the number of runs.
#' Each row in each matrix corresponds to a country and each column to a time.
#' the rownames of the matrices are the country names.
#' @param country_name name of country for which to plot time series
#' @param y_label string. y-axis label
#' @param thin_every an integer.  If 0, ignore.  Otherwise, thin each matrix
#' every this many columns
#' @param thin_integer a logical. If TRUE, thin each matrix to integer values of 
#' times.  If FALSE, ignore
#' @param thin_by_sum logical.  IF TRUE, thin by adding together columns (e.g.
#' if thin_every = 2, add togeher columns 1-2, columns 3-4 etc.)  This is useful
#' for thinning quantities such as incidence.  If FALSE, ignore.
#' @return ggplot object with 95% CI and ten runs
#' @import ggplot2
#' @export
plot_country_time_series <- function(output, country_name, y_label, 
                                     thin_every = 0, 
                                     thin_integer = FALSE,
                                     thin_by_sum = FALSE) {
  stopifnot(is_integer_like(thin_every) && thin_every >= 0)
  country_time_series <- lapply(output, function(x) x[country_name,])
  country_time_series <- do.call(rbind, country_time_series)
  rownames(country_time_series) <- seq_len(nrow(country_time_series))
  
  if(thin_every > 0 || thin_integer) {
    country_time_series <- thin_time_series(country_time_series, 
                                            thin_every = thin_every, 
                                            thin_integer = thin_integer,
                                            thin_by_sum = thin_by_sum)
  }
  
  quantiles <- c(.025, .975)
  quantile_df <- t(apply(country_time_series, 2, function(x) quantile(x, quantiles)))
  quantile_df <- as.data.frame(quantile_df)
  quantile_df$time <- as.numeric(rownames(quantile_df))
  colnames(quantile_df) <- c("lower", "upper", "time")
  
  max_plot_runs <- 10
  
  if(nrow(country_time_series) > max_plot_runs) {
    country_time_series <- country_time_series[seq_len(max_plot_runs), ]
  }

  country_time_series <- melt(country_time_series)
  colnames(country_time_series) <- c("run", "time", "y")
  country_time_series$run <- factor(country_time_series$run)
  
  ggplot(country_time_series) + 
    geom_line(aes(x = time, y = y, group = run)) + 
    geom_ribbon(data = quantile_df, 
                mapping = aes(x = time, ymin = lower, ymax = upper), 
                fill = "red", 
                alpha = .2) +
    scale_x_continuous("Time (days)", expand = c(0, 0)) +
    scale_y_log10(y_label, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max(country_time_series$time)),
                    ylim = c(1, max(country_time_series$y) * 1.1)) +
    theme_bw() +
    theme(text = element_text(size = 24))

}

#' make a string to label the y-axis
#' 
#' @param y_label_base string
#' @param thin_every an integer.  If 0, ignore.  Otherwise, the matrix was thinned
#' every this many columns
#' @param thin_integer a logical. If TRUE, the matrix was thinned to integer values of 
#' times.  If FALSE, ignore
#' @return a string
make_y_label <- function(y_label_base, thin_every, thin_integer) {
  if(thin_every %in% c(0, 1) || thin_integer) {
    y_label <- paste0("Daily ", y_label_base)
  } else if(thin_every == 7) {
    y_label <- paste0("Weekly ", y_label_base)
  } else {
    y_label <- paste0(y_label_base, " every ", thin_every, " days")
  }
  y_label
}
