#' tssDist
#'
#' Plot distance to TSS distribution
#' @param tss_pos File containing "gene chrom tss" information [Default 'data/tss_positions.txt']
#' @import ggplot2
#' @keywords tss
#' @export

tssDist <- function(tss_pos="data/tss_positions.txt"){
  tss_locations<-read.delim(tss_pos, header = T)
  data<-getData()
  
  fun <- function(p) {
    
    index<-which.min(abs(tss_locations$tss - p))
    closestTss<-tss_locations$tss[[index]]
    dist<-(closestTss-p)
    list(closest=closestTss, distance2nearest=dist)
  }
  
  dist2tss<-lapply(data$pos, fun)
  
  dist2tss<-do.call(rbind, dist2tss)
  dist2tss<-as.data.frame(dist2tss)
  dist2tss<-arrange(dist2tss,(as.numeric(distance2nearest)))
  
  dist2tss$distance2nearest<-as.numeric(dist2tss$distance2nearest)
  
  p<-ggplot(dist2tss)
  p<-p + geom_density(aes(distance2nearest), fill = "lightblue", alpha = 0.6)
  p<-p + scale_x_continuous("Distance to TSS", limits=c(-10000, 10000))
  p
  
  # p<-ggplot(dist2tss)
  # p<-p + geom_histogram(aes(distance2nearest), fill = "lightblue", alpha = 0.6, bins=500)
  # p<-p + scale_x_continuous("Distance to TSS", limits=c(-10000, 10000))
  # p
  
  #return(dist2tss)
}