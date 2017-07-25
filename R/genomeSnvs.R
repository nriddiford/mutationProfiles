#' genomeSnvs
#'
#' Plot snvs accross genome
#' @import ggplot2
#' @keywords genome
#' @export


genomeSnvs <- function(){
  data<-getData()
  data<-filter(data, chrom != "Y" & chrom != "4")
  p<-ggplot(data)
  p<-p + geom_point(aes(pos/1000000, sample, colour = trans))
  # p<-p + guides(color = FALSE)
  p<-p + theme(axis.text.x = element_text(angle=45, hjust = 1))
  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  
  p
}


