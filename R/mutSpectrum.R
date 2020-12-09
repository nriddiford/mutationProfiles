#' mutSpectrum
#'
#' Plots the mutations spectrum for all samples combined
#' @import ggplot2
#' @keywords spectrum
#' @export
mutSpectrum <- function(..., snv_data=NULL, write=FALSE, max_y=25, collapse_tris = FALSE){
  if(missing(snv_data)){
    snv_data<-getData(...)
  }

  snv_data <- snv_data %>%
    dplyr::filter(...)

  # Should correct for genome wide occurances...
  # comnined_data <- snv_data %>%
  #   dplyr::group_by(decomposed_tri, grouped_trans) %>%
  #   dplyr::tally()
  #
  #
  # scaling_factor <-triFreq()
  # plyr::join(comnined_data, )


  cat("Showing global contribution of tri class to mutation load", "\n")

  p <- ggplot(snv_data)
  if(collapse_tris){
    p <- p + geom_bar(aes(x = grouped_trans, y = (..count..)/sum(..count..), group = grouped_trans, fill = grouped_trans), alpha=0.7, position="dodge",stat="count")
  } else {
    p <- p + geom_bar(aes(x = decomposed_tri, y = (..count..)/sum(..count..), group = decomposed_tri, fill = grouped_trans), alpha=0.7, position="dodge",stat="count")
  }
  p <- p + scale_y_continuous("Contribution to mutation load", limits = c(0, max_y/100), breaks=seq(0,max_y/100,by=0.025), labels=paste0(seq(0,max_y,by=2.5), "%"), expand = c(0.0, .0005))
  p <- p + scale_x_discrete("Genomic context", expand = c(.005, .005))
  p <- p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5, size=9),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=15),
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

