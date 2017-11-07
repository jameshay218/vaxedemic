#' Plot total arrivals and departures
#' 
#' Takes a melted data frame and then plots the aggregated number of flights/connections
#' between each pair of the Destination and Origin column
#' @param data a data frame which MUST have a column for Destination, Origin and Volume
#' @param uselog if TRUE, plots the log value
#' @return a ggplot2 object
#' @export
plot_total_journies <- function(data, uselog=FALSE){

  total_arrivals_to <- plyr::ddply(data[,c("Destination","Volume")],~Destination,function(x) sum(x$Volume))
  total_arrivals_from <- plyr::ddply(data[,c("Origin","Volume")],~Origin,function(x) sum(x$Volume))
  
  if(uselog){
    total_arrivals_to$V1 <- log(total_arrivals_to$V1)
    total_arrivals_from$V1 <- log(total_arrivals_from$V1)
  }
  
  p_arrivals_to <- ggplot(total_arrivals_to) + 
    geom_bar(aes(x=Destination,y=V1),stat="identity")+
    theme_bw()+
    ylab("Number of arrivals") +
    xlab("Destination country") +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
          axis.text.y=element_text(size=6),
          panel.grid = element_blank()) +
    ggtitle("Total number of arrivals coming *to* country")
  
  p_arrivals_from <- ggplot(total_arrivals_from) + 
    geom_bar(aes(x=Origin,y=V1),stat="identity")+
    theme_bw() +
    ylab("Number of departures") +
    xlab("Source country") +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
          axis.text.y=element_text(size=6),
          panel.grid = element_blank()) +
    ggtitle("Total number of arrivals coming *from* country")
  
  p_arrivals <- cowplot::plot_grid(p_arrivals_to,p_arrivals_from,ncol=1)
  p_arrivals
}

#' Log connectivity heatmap
#' 
#' Plots a heatmap of log connectivity between two countries
#' @param data a data frame which MUST have a column for Destination, Origin and logVolume
#' @return a ggplot2 object
#' @export
plot_heatmap_travel <- function(data){
  p1 <- ggplot(data) + 
    geom_bin2d(aes(x=Destination,y=Origin, fill=logVolume),stat="identity") +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=4),
          axis.text.y=element_text(size=4))
  p1
}