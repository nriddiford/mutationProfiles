#' geneHit
#'
#' Show top hit genes
#' @import dplyr
#' @keywords gene
#' @export


geneHit <- function(){
  data<-getData()
  
  data<-filter(data, gene != "intergenic")
  
  hit_count<-sort(table(unlist(data$gene)), decreasing = T)
  head(hit_count)

}