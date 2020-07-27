#' getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords get
#' @import dplyr plyr
#' @export
#' @return Dataframe
getData <- function(..., infile = system.file("extdata", "annotated_snvs.txt", package="mutationProfiles"),
                    expression_source = 'flygut',
                    expression_data='/Users/Nick_curie/Documents/Curie/Data/RNA-Seq_data/Buchon_summary_ISCs.txt',
                    type='snv',
                    attach_info='data/samples_names_conversion.txt'){

  cat("Filters applied:\n")
  input_list <- as.list(substitute(list(...)))
  lapply(X=input_list, function(x) {str(x);summary(x)})
  snv_data<-read.delim(infile, header = T)
  cat("Reading SNVs from ", infile, "\n")

  if(type=='snv'){
    colnames(snv_data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "af", "caller", "variant_type", "status", "snpEff_anno", "feature", "gene", "id")
    # snv_data$dups<-duplicated(snv_data[,1:3])
    # snv_data<- dplyr::mutate(snv_data, caller = ifelse(dups == "TRUE", 'varscan2_mutect2' , as.character(caller)))
  }

  if(type=='indel'){
    colnames(snv_data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "type", "decomposed_tri", "af", "caller", "variant_type", "status", "snpEff_anno", "feature", "gene", "id")
    snv_data$alt <- gsub('\\+|-', '', snv_data$alt)
    snv_data$alt <- gsub("\\'", '', snv_data$alt)
  }

  if(file.exists(attach_info)){
    name_conversion <- read.delim(attach_info, header=F)
    cat("Attaching assay information to data\n")
    colnames(name_conversion) <- c("sample", "sample_short", "sample_paper", "sex", "assay")
    snv_data <- plyr::join(snv_data, name_conversion, "sample", type = 'left') %>%
      # dplyr::rename(sample_old = sample,
      #               sample = sample_paper) %>%
      dplyr::select(sample_paper, everything())
  }
  if(expression_source == 'flygut'){
    cat("Reading expression data from source: 'flygut [Buchon]'\n\n")
    expression_data = read.delim('/Users/Nick_curie/Documents/Curie/Data/RNA-Seq_data/Buchon_summary_ISCs.txt')
    colnames(expression_data) <- c('id', 'symbol', 'name', 'isc', 'eb', 'ec', 'ee', 'vm')
    seq_data <- expression_data %>%
      dplyr::mutate(fpkm = isc) %>%
      dplyr::select(id, fpkm)
  } else{
    cat("Reading expression data from source: 'Dutta'\n\n")
    expression_data = system.file("extdata", "isc_genes_rnaSeq.txt")
    # Read in tissue specific expression data
    seq_data<-read.csv(header = F, expression_data)
    colnames(seq_data)<-c('id', 'fpkm')
  }

  snv_data <- plyr::join(snv_data,seq_data, "id", type = 'left')

  snv_data <- snv_data %>%
    dplyr::filter(...) %>%
    dplyr::mutate(fpkm = ifelse(is.na(fpkm), 0, round(fpkm, 1))) %>%
    dplyr::mutate(pos = as.numeric(pos),
                  af = as.double(af)) %>%
    dplyr::mutate(cell_fraction = ifelse(chrom %in% c('X', 'Y'), af,
                                         ifelse(af*2>1, 1, af*2))) %>%
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
  # print(levels(snv_data$sample))
  dput(as.character(levels(snv_data$sample)))
}
