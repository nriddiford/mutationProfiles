#' notchSnvs
#'
#' Plot snvs around Notch by sample
#' @import ggplot2
#' @keywords notch
#' @export

notchSnvs <- function(){
  data<-getData()
  data<-filter(data, chrom == "X", pos >= 3000000, pos <= 3300000)
  
  if(nrow(data) == 0){
    stop("There are no snvs in Notch. Exiting\n")
  }
  
  p<-ggplot(data)
  p<-p + geom_point(aes(pos/1000000, sample, colour = trans, size = 2))
  p<-p + guides(size = FALSE, sample = FALSE)
  p<-p + cleanTheme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(3,3.3,by=0.05), limits=c(3, 3.301))
  p<-p + annotate("rect", xmin=3.000000, xmax=3.134532, ymin=0, ymax=0.1, alpha=.2, fill="green")
  p<-p + annotate("rect", xmin=3.134870, xmax=3.172221, ymin=0, ymax=0.1, alpha=.2, fill="skyblue")
  p<-p + annotate("rect", xmin=3.176440, xmax=3.300000, ymin=0, ymax=0.1, alpha=.2, fill="red")
  
  p
}