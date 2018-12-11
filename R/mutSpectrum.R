#' mutSpectrum
#'
#' Plots the mutations spectrum for all samples combined
#' @import ggplot2
#' @keywords spectrum
#' @export

mutSpectrum <- function(){
  data<-getData()
  cat("Showing global contribution of tri class to mutation load", "\n")
 
  p<-ggplot(data)
  p<-p + geom_bar(aes(x = decomposed_tri, y = (..count..)/sum(..count..), group = decomposed_tri, fill = grouped_trans), position="dodge",stat="count")
  p<-p + scale_y_continuous("Relative contribution to mutation load", expand = c(0.0, .0005))
  p<-p + scale_x_discrete("Genomic context", expand = c(.005, .005))
  p<-p + cleanTheme() + 
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15)
          )
  p<-p + labs(fill="Mutation class")
  p<-p + facet_wrap(~grouped_trans, ncol = 3, scale = "free_x" )
  
  mut_spectrum<-paste("mutation_spectrum.pdf")
  cat("Writing file", mut_spectrum, "\n")
  ggsave(paste("plots/", mut_spectrum, sep=""), width = 20, height = 10)
  p
}