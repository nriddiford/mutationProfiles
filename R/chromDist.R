#' chromDist
#'
#' Plot genome-wide snv distribution
#' @import ggplot2
#' @keywords distribution
#' @export


chromDist <- function(..., snv_data=NULL, object=NA, notch=0, write=FALSE){
  if(missing(snv_data)){
    snv_data<-getData(...)
  }
  ext<-'.pdf'
  if(is.na(object)){
    object<-'grouped_trans'
    cols<-setCols(snv_data, "grouped_trans")
  }

  if(notch){
    snv_data<-exclude_notch(...)
    ext<-'_excl.N.pdf'
  }

  cat("Plotting snvs by", object, "\n")

  p<-ggplot(snv_data)
  p<-p + geom_histogram(aes(pos/1000000, fill = get(object)), binwidth=0.1, alpha = 0.8)
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33),expand = c(0.01, 0.01))
  p<-p + scale_y_continuous("Number of snvs", expand = c(0.01, 0.01))
  p<-p + cleanTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15)
    )

  if (object == 'grouped_trans'){
    p<-p + cols
  }
  if(write){
    chrom_outfile<-paste("snv_dist_genome_by_", object, ext, sep = "")
    cat("Writing file", chrom_outfile, "\n")
    ggsave(paste("plots/", chrom_outfile, sep=""), width = 20, height = 10)
  }
  p
}
