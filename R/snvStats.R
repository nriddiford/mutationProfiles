#' snvStats
#'
#' Calculate some basic stats for snv snv_data
#' @import dplyr
#' @keywords stats
#' @export

snvStats <- function(){
  snv_data<-getData()
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
