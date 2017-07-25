#' genTris
#'
#' This function returns all possible trinucleotide combinations
#' @keywords trinucleotides
#' @export
#' @return Character string containing all 96 trinucleotides
#' genTris()


genTris <- function(){
  all.tri = c()
    for(i in c("A", "C", "G", "T")){
      for(j in c("C", "T")){
        for(k in c("A", "C", "G", "T")){
          if(j != k){
            for(l in c("A", "C", "G", "T")){
              tmp = paste(i, "[", j, ">", k, "]", l, sep = "")
              all.tri = c(all.tri, tmp)
            }
          }
        }
      }
    }
  all.tri <- all.tri[order(substr(all.tri, 3, 5))]
  return(all.tri)
}