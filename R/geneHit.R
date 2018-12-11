#' geneHit
#'
#' Show top hit genes
#' @import dplyr
#' @keywords gene
#' @param n Show top n hits [Default 10]
#' @export
geneHit <- function(..., n=10){
  snv_data<-getData(...)
  snv_data<-filter(snv_data, gene != "intergenic")

  hit_count<-as.data.frame(sort(table(unlist(snv_data$gene)), decreasing = T))

  colnames(hit_count)<- c("gene", "count")
  head(hit_count, n)
}
