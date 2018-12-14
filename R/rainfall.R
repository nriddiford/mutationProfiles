#' rainfall
#'
#' Plot log10 distances between snvs as rainfall plot
#' @import ggplot2
#' @keywords rainfall
#' @export

rainfall <- function(..., snv_data=NULL){
  if(missing(snv_data)){
    snv_data<-getData(...)
  }

  distances <- do.call(rbind, lapply(split(snv_data[order(snv_data$chrom, snv_data$pos),], snv_data$chrom[order(snv_data$chrom, snv_data$pos)]),
                                     function(a)
                                       data.frame(a,
                                                  dist=c(diff(a$pos), NA),
                                                  logdist = c(log10(diff(a$pos)), NA))
                                     )
                       )

  distances$logdist[is.infinite(distances$logdist)] <- 0
  distances<-filter(distances, chrom != 4)

  p <- ggplot(distances)
  p <- p + geom_point(aes(pos/1000000, logdist, colour = grouped_trans))
  p <- p + cleanTheme() +
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          strip.text = element_text(size=20)
    )

  p <- p + facet_wrap(~chrom, scale = "free_x", ncol = 6)
  #p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  p <- p + scale_x_continuous("Mbs", breaks = seq(0,max(distances$pos),by=10))

  rainfall_out<-paste("rainfall.pdf")
  cat("Writing file", rainfall_out, "\n")
  ggsave(paste("plots/", rainfall_out, sep=""), width = 20, height = 5)

  p
}
