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
