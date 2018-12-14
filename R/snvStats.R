#' snvStats
#'
#' Calculate some basic stats for snv snv_data
#' @import dplyr
#' @keywords stats
#' @export

snvStats <- function(..., snv_data=NULL){
  if(missing(snv_data)){
    snv_data<-getData(...)
  }

  cat("sample", "snvs", sep='\t', "\n")
  rank<-sort(table(snv_data$sample), decreasing = TRUE)
  rank<-as.array(rank)

  total=0

  scores=list()
  for (i in 1:nrow(rank)){
    cat(names(rank[i]), rank[i], sep='\t', "\n")
    total<-total + rank[i]
    scores[i]<-rank[i]
  }


  cat('--------------', '\n')
  scores<-unlist(scores)

  mean<-as.integer(mean(scores))
  med<-as.integer(median(scores))

  cat('total', total, sep='\t', '\n')
  cat('samples', nrow(rank), sep='\t', '\n')

  cat('--------------', '\n')
  cat('mean', mean, sep='\t', '\n')
  cat('median', med, sep='\t', '\n')

  cat('\n')
  all_ts<-nrow(filter(snv_data, trans == "A>G" | trans == "C>T" | trans == "G>A" | trans == "T>C"))
  all_tv<-nrow(filter(snv_data, trans != "A>G" & trans != "C>T" & trans != "G>A" & trans != "T>C"))
  ts_tv<-round((all_ts/all_tv), digits=2)
  cat("ts/tv = ", ts_tv,  sep='', '\n')

}

#' splot_snvs
#'
#' Plot snv count per sample
#' @import dplyr
#' @keywords count
#' @export
plot_snvs <- function(..., snv_data=NULL){
  if(missing(snv_data)){
    snv_data<-getData(...)
  }

  calls_by_sample <- snv_data %>%
    dplyr::group_by(sample) %>%
    dplyr::tally()

  calls_by_sample <- transform(calls_by_sample, sample = reorder(sample, -n))

  p <- ggplot(calls_by_sample, aes(sample, n))
  p <- p + geom_bar(stat='identity')
  p <- p + scale_y_continuous("Number of calls", limits=c(0,max(calls_by_sample$n)), breaks=seq(0,max(calls_by_sample$n), by=100))
  p <- p + scale_x_discrete("Sample")
  p <- p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
          axis.text = element_text(size=20), axis.title = element_text(size=20)
    )

  p


}
