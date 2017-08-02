#' featuresHit
#'
#' Show top hit features
#' @import ggplot2
#' @keywords features
#' @export


featuresHit <- function(){
  data<-getData()
  
  # To condense exon counts into "exon"
  data$feature<-as.factor(gsub("exon_.*", "exon", data$feature))
  
  # Reoders descending
  data$feature<-factor(data$feature, levels = names(sort(table(data$feature), decreasing = TRUE)))
  
  #cols<-set_cols(data, "feature")
  
  p<-ggplot(data)
  p<-p + geom_bar(aes(feature, fill = feature))
  #p<-p + cols
  p<-p + cleanTheme() +
    theme(axis.title.x=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous(expand = c(0.01, 0.01))
  
  features_outfile<-paste("hit_features_count.pdf")
  cat("Writing file", features_outfile, "\n")
  ggsave(paste("plots/", features_outfile, sep=""), width = 20, height = 10)
  
  p
}