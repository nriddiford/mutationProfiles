#' chromDist
#'
#' Plot genome-wide snv distribution
#' @import ggplot2
#' @keywords distribution
#' @param snv_data A dataframe of snvs [Default: NULL]
#' @param notch Filter in/out Notch region [Default: FALSE]
#' @param write Write plot to file 'plots/snv_dist_genome_by.png'? [Default: FALSE]
#' @param binsize Control the density adjustment if density=TRUE, or the binwidth if density=FALSE. For histograms a value of 0.1 corresponds to 1/10th of a Mb [Default:0.1]
#' @param density Plot as density? [Default: TRUE]
#' @export
chromDist <- function(..., snv_data=NULL, notch=FALSE, write=FALSE, binsize=0.1, density=TRUE){
  if(notch){
    snv_data<-exclude_notch(...)
    ext<-'_excl.N.png'
  } else if(missing(snv_data)){
    snv_data<-getData(...)
  }

  cols<-setCols(snv_data, "grouped_trans")
  blueBar <- '#3B8FC7'

  cat("Plotting snv distribution across genome", "\n")

  p <- ggplot(snv_data)
  if(density){
    p <- p + geom_density(aes(pos/1000000, fill=chrom), alpha = 0.6, adjust=binsize)
    p <- p + scale_y_continuous("Density", expand = c(0.01, 0.01))

  } else{
    p <- p + geom_histogram(aes(pos/1000000, fill=blueBar), binwidth=binsize, alpha = 0.8)
    p <- p + scale_y_continuous("Number of SNVs", expand = c(0.01, 0.01))
    p <- p + scale_fill_identity()
  }
  # p<-p + geom_histogram(aes(pos/1000000, fill = grouped_trans), binwidth=0.1, alpha = 0.8)
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Genomic position (Mbs)", expand = c(0.01, 0.01))
  p<-p + cleanTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15)
    )
  # p <- p + cols
  if(write){
    chrom_outfile<-"snv_dist_genome_by.png"
    cat("Writing file", paste0("plots/", chrom_outfile, "\n"))
    ggsave(paste0("plots/", chrom_outfile), width = 20, height = 10)
  }
  p
}
