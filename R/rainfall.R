#' rainfall
#'
#' Plot log10 distances between snvs as rainfall plot
#' @import ggplot2
#' @keywords rainfall
#' @export

rainfall <- function(){
  data<-getData()

  distances<-do.call(rbind, lapply(split(data[order(data$chrom, data$pos),], data$chrom[order(data$chrom, data$pos)]),
                        function(a) 
                          data.frame(a,
                                     dist=c(diff(a$pos), NA),
                                     logdist = c(log10(diff(a$pos)), NA))
                        )
                     )

  distances$logdist[is.infinite(distances$logdist)] <- 0
  distances<-filter(distances, chrom != 4)

  p<-ggplot(distances)
  p<-p + geom_point(aes(pos/1000000, logdist, colour = grouped_trans))
  # p<-p + guides(color = FALSE)
  p<-p + theme(axis.text.x = element_text(angle=45, hjust = 1))
  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  
  p
}