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

