#' setCols
#'
#' Show top hit genes
#' @import RColorBrewer
#' @param df Dataframe [Required]
#' @param col Column of dataframe. Colours will be set to levels(df$cols) [Required]
#' @keywords cols
#' @export


setCols <- function(df, col){
  names<-levels(df[[col]])
  cat("Setting colour levles:", names, "\n")
  level_number<-length(names)
  mycols<-brewer.pal(level_number, "Set2")
  names(mycols) <- names
  colScale <- scale_fill_manual(name = col,values = mycols)
  return(colScale)
}
