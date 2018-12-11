#' mutSpectrum
#'
#' Plots the mutations spectrum for all samples combined
#' @import ggplot2
#' @keywords spectrum
#' @export
mutSpectrum <- function(write=FALSE, max_y=25){
  snv_data<-getData()
  cat("Showing global contribution of tri class to mutation load", "\n")

  p <- ggplot(snv_data)
  p <- p + geom_bar(aes(x = decomposed_tri, y = (..count..)/sum(..count..), group = decomposed_tri, fill = grouped_trans), position="dodge",stat="count")
  p <- p + scale_y_continuous("Contribution to mutation load", limits = c(0, max_y/100), breaks=seq(0,max_y/100,by=0.025), labels=paste0(seq(0,max_y,by=2.5), "%"), expand = c(0.0, .0005))
  p <- p + scale_x_discrete("Genomic context", expand = c(.005, .005))
  p <- p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 90, hjust=1),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15)
    )
  p <- p + facet_wrap(~grouped_trans, ncol = 6, scale = "free_x" )
  p <- p + guides(grouped_trans = FALSE)
  if(write){
    mut_spectrum<-paste("mutation_spectrum.pdf")
    cat("Writing file", mut_spectrum, "\n")
    ggsave(paste("plots/", mut_spectrum, sep=""), width = 20, height = 5)
  }
  p
}

