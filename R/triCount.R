#' triFreq
#'
#' This function counts the number of times each triunucleotide is found in a supplied genome
#' @param genome BS.genome file defaults to BSgenome.Dmelanogaster.UCSC.dm6
#' @param count Output total counts instead of frequency if set [Default no]
#' @import dplyr
#' @keywords trinucleotides
#' @export
#' @return Dataframe of trinucs and freqs (or counts if count=1)
triFreq <- function(genome=NA, count=NA){
  if(is.na(genome)){
    cat("No genome specfied, defaulting to 'BSgenome.Dmelanogaster.UCSC.dm6'\n")
    library(BSgenome.Dmelanogaster.UCSC.dm6, quietly = TRUE)
    genome <- BSgenome.Dmelanogaster.UCSC.dm6
  }

  params <- new("BSParams", X = Dmelanogaster, FUN = trinucleotideFrequency, exclude = c("M", "_"), simplify = TRUE)
  data<-as.data.frame(bsapply(params))
  data$genome<-as.integer(rowSums(data))
  data$genome_adj<-(data$genome*2)

  if(!is.na(count)){
    tri_count<-data['genome_adj']
    tri_count<-cbind(tri = rownames(tri_count), tri_count)
    colnames(tri_count) <- c("tri", "count")
    rownames(tri_count) <- NULL
    return(tri_count)
  }
  else{
    data$x <- (1/data$genome)
    scaling_factor<-data['x']
    return(scaling_factor)
  }

}

#' A function to calculate trinucleotide frequences from a supplied fasta file
#' @param fastaFile A fasta file (genome, exome)
#' @import dplyr Biostrings
#' @export
#' @return Dataframe of trinucs and freqs (or counts if count=1)
getTriFromFasta <- function(fastaFile = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta'){
  x <- readDNAStringSet(fastaFile)
  triCount <- trinucleotideFrequency(x, simplify.as = "collapsed")
  d <- as.data.frame(triCount)
  d$genome_adj<-(d$tri*2)
  d$x <- 1/d$tri

  data_frame(
    data = triCount,
    cat = names(triCount))

  d <- triCount %>%
    data.frame() %>%
    mutate(genome_adj = tri *2) %>%
    mutate(x = 1/tri) %>%
    select(x)

  row.names(d) <- names(triCount)
  return(d)
}
