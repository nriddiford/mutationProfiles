#' featuresHit
#'
#' Show top hit features
#' @import ggplot2
#' @keywords features
#' @export

featuresHit <- function(..., snv_data=NULL, write=FALSE){
  if(missing(snv_data)){
    snv_data<-getData(...)
  }

  # To condense exon counts into "exon"
  snv_data$feature<-as.factor(gsub("exon_.*", "exon", snv_data$feature))

  # Reoders descending
  snv_data$feature<-factor(snv_data$feature, levels = names(sort(table(snv_data$feature), decreasing = TRUE)))

  snv_data <- snv_data %>%
    dplyr::group_by(feature) %>%
    dplyr::add_tally() %>%
    ungroup() %>%
    dplyr::filter(n >= 5) %>%
    droplevels()
  #cols<-setCols(snv_data, "feature")

  p <- ggplot(snv_data)
  p <- p + geom_bar(aes(feature, fill = feature))
  #p<-p + cols
  p <- p + cleanTheme() +
    theme(axis.title.x=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p <- p + scale_x_discrete(expand = c(0.01, 0.01))
  p <- p + scale_y_continuous(expand = c(0.01, 0.01))

  # colour to a pub palette:
  # p<-p + ggpar(p, palette = 'jco')

  if(write){
    features_outfile<-paste("hit_features_count.pdf")
    cat("Writing file", features_outfile, "\n")
    ggsave(paste("plots/", features_outfile, sep=""), width = 20, height = 10)
  }
  p
}
