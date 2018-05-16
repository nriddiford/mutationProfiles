#' getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords get
#' @import dplyr
#' @export
#' @return Dataframe

getData <- function(infile = "data/annotated_snvs.txt"){
  data<-read.delim(infile, header = F)

  colnames(data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "caller", "type", "feature", "gene", "")
  levels(data$type) <- tolower(levels(data$type))
  #data <- filter(data, type == 'germline')

  #filter on chroms
  # data<-filter(data, chrom != "Y" & chrom != "4")
  #filter out samples
  # data<-filter(data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  data<-droplevels(data)
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(data)
}

