#' rainfall
#'
#' Plot log10 distances between snvs as rainfall plot
#' @import ggplot2 dplyr
#' @keywords rainfall
#' @export

rainfall <- function(..., snv_data=NULL, write=FALSE, title=NULL, exclude_chroms=c("4", "Y")){
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
  distances <- dplyr::filter(distances, !chrom %in% exclude_chroms)


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
  if(!missing(title)){
    p <- p + ggtitle(title)
  }

  if(write){
    rainfall_out<-paste0(title, "_rainfall.png")
    cat("Writing file", rainfall_out, "\n")
    ggsave(paste("plots/", rainfall_out, sep=""), width = 20, height = 5)
  }

  p

}
