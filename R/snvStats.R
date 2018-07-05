#' snvStats
#'
#' Calculate some basic stats for snv data
#' @import dplyr
#' @keywords stats
#' @export
snvStats <- function(){
  data<-getData()
  cat("sample", "snvs", sep='\t', "\n")
  rank<-sort(table(data$sample), decreasing = TRUE)
  rank<-as.array(rank)

  for (i in 1:nrow(rank)){
    cat(names(rank[i]), rank[i], sep='\t', "\n")
  }

  all_ts<-nrow(filter(data, trans == "A>G" | trans == "C>T" | trans == "G>A" | trans == "T>C"))
  all_tv<-nrow(filter(data, trans != "A>G" & trans != "C>T" & trans != "G>A" & trans != "T>C"))
  ts_tv<-all_ts/all_tv
  cat("ts/tv =", ts_tv)
}
