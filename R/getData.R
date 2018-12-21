#' getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords get
#' @import dplyr plyr
#' @export
#' @return Dataframe
getData <- function(..., infile = "data/annotated_snvs.txt", exclude=TRUE, expression_data='data/isc_genes_rnaSeq.csv', type='snv'){

  cat("Filters applied:\n")
  input_list <- as.list(substitute(list(...)))
  lapply(X=input_list, function(x) {str(x);summary(x)})

  snv_data<-read.delim(infile, header = T)
  if(type=='snv'){
    colnames(snv_data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "af", "caller", "variant_type", "status", "snpEff_anno", "feature", "gene", "id")
    snv_data$dups<-duplicated(snv_data[,1:3])
    snv_data<- dplyr::mutate(snv_data, caller = ifelse(dups == "TRUE", 'varscan2_mutect2' , as.character(caller)))
  }
  if(type=='indel'){
    colnames(snv_data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "type", "decomposed_tri", "af", "caller", "variant_type", "status", "snpEff_anno", "feature", "gene", "id")
    snv_data$alt <- gsub('\\+|-', '', snv_data$alt)
    snv_data$alt <- gsub("\\'", '', snv_data$alt)
  }
  # Read in tissue specific expression data
  seq_data<-read.csv(header = F, expression_data)
  colnames(seq_data)<-c('id', 'fpkm')

  snv_data <- plyr::join(snv_data,seq_data,"id", type = 'left')

  excluded_samples <- c()
  if(exclude){
    excluded_samples <- c("A373R7", "A512R17", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  }
  snv_data <- snv_data %>%
    dplyr::filter(!sample %in% excluded_samples) %>%
    dplyr::mutate(fpkm = ifelse(is.na(fpkm), 0, round(fpkm, 1))) %>%
    dplyr::mutate(pos = as.numeric(pos),
                  af = as.double(af)) %>%
    dplyr::mutate(cell_fraction = ifelse(chrom %in% c('X', 'Y'), af,
                                         ifelse(af*2>1, 1, af*2))) %>%
    dplyr::filter(...) %>%
    droplevels()

  return(snv_data)
}


#' showSamples
#'
#' A helper function to print the sample names
#' @param infile File to process [Required]
#' @keywords samples
#' @import dplyr
#' @export

showSamples <- function(infile = "data/annotated_snvs.txt"){
  snv_data<-read.delim(infile, header = T)
  print(levels(snv_data$sample))
}
