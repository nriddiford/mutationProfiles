#' getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords get
#' @import dplyr plyr
#' @export
#' @return Dataframe
getData <- function(infile = "data/annotated_snvs.txt", expression_data='data/isc_genes_rnaSeq.csv'){
  snv_data<-read.delim(infile, header = F)
  colnames(snv_data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "a_freq", "caller", "variant_type", "status", "snpEff_anno", "feature", "gene", "id")

  # Read in tissue specific expression data
  seq_data<-read.csv(header = F, expression_data)
  colnames(seq_data)<-c('id', 'fpkm')

  snv_data <- plyr::join(snv_data,seq_data,"id", type = 'left')

  snv_data$fpkm <- ifelse(is.na(snv_data$fpkm), 0, round(as.numeric(snv_data$fpkm), 1))

  # Order by FPKM
  snv_data<- dplyr::arrange(snv_data, desc(fpkm))

  # Find vars called by both Mu and Var
  # Must also filter one of these calls out...
  snv_data$dups<-duplicated(snv_data[,1:3])
  # snv_data<- dplyr::mutate(snv_data, caller = ifelse(dups == "TRUE", 'varscan2_mutect2' , as.character(caller)))

  ##############
  ## Filters ###
  ##############

  # Filter for calls made by both V and M
  # snv_data<-filter(snv_data, caller == 'mutect2' | caller == 'varscan2_mutect2')

  # Filter for old/new data
  # cat("Filtering for old/new data\n")
  # snv_data <- filter(snv_data, !grepl("^A|H", sample))

  # Filter for genes expressed in RNA-Seq data
  # cat("Filtering out non-expressed genes\n")
  # snv_data<-filter(snv_data, !is.na(fpkm) & fpkm > 0.1)

  # Filter for genes NOT expressed in RNA-Seq data
  # cat("Filtering out expressed genes\n")
  # snv_data<-filter(snv_data, fpkm == 0)

  # Filter on allele freq
  # cat("Filtering on allele frequency\n")
  #snv_data<-filter(snv_data, is.na(a_freq))
  # snv_data<-filter(snv_data, a_freq >= 0.20)

  # Filter out samples
  # snv_data<-filter(snv_data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  # snv_data <- dplyr::filter(snv_data, !sample %in% c("A373R1", "A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9"))
  # snv_data<-filter(snv_data, sample != "A373R11" & sample != 'A373R13')

  # snv_data <- snv_data %>%
  #   filter(sample %in% c("D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2")) %>%
  #   droplevels()
  dir.create(file.path("plots"), showWarnings = FALSE)

  return(snv_data)
}
