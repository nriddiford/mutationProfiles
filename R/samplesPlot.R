#' samplesPlot
#'
#' Plot the snv distribution for each sample
#' @import ggplot2
#' @param count Output total counts instead of frequency if set [Default no]
#' @keywords spectrum
#' @export

samplesPlot <- function(..., snv_data=NULL, count=FALSE, write=FALSE, tricounts=NULL){
  if(missing(snv_data)){
    snv_data<-getData(...)
  }
  mut_class<-c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

  if(missing(tricounts)){
    cat("calculating trinucleotide frequencies in genome\n")
    tricounts <- triFreq()
  }

  frequencies <- snv_data %>%
    dplyr::group_by(sample, grouped_trans) %>%
    dplyr::tally()

  p <- ggplot(snv_data)

  if(!count){
    p<-p + geom_bar(aes(x = grouped_trans, y = (..count..)/sum(..count..), group = sample, fill = sample), position="dodge",stat="count")
    p<-p + scale_y_continuous("Relative contribution to total mutation load", expand = c(0.0, .001))
    tag='_freq'
  }
  else{
    p<-p + geom_bar(aes(x = grouped_trans, y = ..count.., group = sample, fill = sample), position="dodge",stat="count")
    p<-p + scale_y_continuous("Count", expand = c(0.0, .001))
    tag='_count'
  }
  p<-p + scale_x_discrete("Mutation class", limits=mut_class)
  p<-p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 10)
          )
  p<-p + facet_wrap(~sample, ncol = 4, scale = "free_x" )

  if(write){
    samples_mut_spect<-paste("mutation_spectrum_samples", tag, ".pdf", sep = '')
    cat("Writing file", samples_mut_spect, "\n")
    ggsave(paste("plots/", samples_mut_spect, sep=""), width = 20, height = 10)
  }
  p
}
