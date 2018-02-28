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
      
    ggplot(I_aggregated[I_aggregated$Location %in% unique(I_aggregated$Location)[n_countries],],aes(x=variable,y=x/X,col=Age)) +
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

#' @export
plot_peak_times <- function(res, labels, regionDat, latitudeDat){
    peakTimes = do.call("cbind",res)
    summaryPeaks <- t(apply(peakTimes, 1, function(x) c(mean(x),quantile(x, c(0.025,0.5,0.975)))))
    colnames(summaryPeaks) <- c("mean","lower95","median","upper95")
    summaryPeaks <- data.frame(Location=unique(labels$Location), summaryPeaks)
    dat <- merge(summaryPeaks,regionDat[,c("Location","region")])
    dat <- merge(dat,latitudeDat)
    tmp <- unique(dat[,c("Location","latitude")])
    my_latitudes <- order(tmp$latitude)
    dat$Location <- factor(dat$Location, levels=tmp$Location[my_latitudes])
    dat$Hemisphere <- ifelse(dat$latitude < 0,"Southern","Northern")
   
    p <- ggplot(dat) +
        geom_errorbarh(aes(y=Location,x=median,xmax=upper95,xmin=lower95, col=Hemisphere)) +
        geom_point(aes(y=Location,x=median),size=0.5) +
        scale_x_continuous(limits=c(0,365)) +
        xlab("Peak time (days), median and 95% quantiles") +
        facet_grid(region~.,scales="free_y", space="free",switch="both") +
        theme(axis.text.y=element_text(size=6)) + theme_bw()
    return(list("plot"=p,"data"=dat))
}
