selctionDetect <- function(..., snv_data = NULL, refDB = "data/dm6_refcds.rda", per_sample = F, onlySample=NULL){
  if(missing(snv_data)) snv_data<-getData(...)

  if(!missing(onlySample)) {
    cat("Running for sample:", onlySample, "\n")
    snv_data <- snv_data %>% dplyr::filter(sample == onlySample) %>% droplevels()
  }

  if(!per_sample) return(rundNdS(snv_data, refDB))

  sample_names <- levels(snv_data$sample)
  outList <- list()
  for(i in 1:length(sample_names)){
    cat("Running per sample:", sample_names[i], "\n")
    # snv_data <- snv_data %>% dplyr::filter(sample == sample_names[i]) %>% droplevels()
    cat(nrow(snv_data[snv_data$sample==sample_names[i],]), "mutations in sample\n")
    tryCatch(
      expr = {
        out <- rundNdS(snv_data[snv_data$sample==sample_names[i],], refDB)
        outList[sample_names[i]] <- out
      },
      error = function(e){
        print(e)
      }
    )
  }
  # return(outList)
  return(compressTable(outList))
}


rundNdS <- function(data, refDB){
  mutations <- data %>%
    dplyr::select(sample, chrom, pos, ref, alt)

  colnames(mutations) <- c('sampleID', 'chr', 'pos', 'ref', 'mut')

  dndsout = dndscv(mutations, refdb = refDB)

  return(dndsout)
}

compressTable <- function(l){
  data <- as.data.frame(do.call(rbind, l))
  data$sample <- rownames(data)
  rownames(data) <- NULL

  data <- data %>%
    dplyr::mutate(sample = gsub("\\..*", "", sample)) %>%
    dplyr::select(sample, name, mle, cilow, cihigh)

  return(data)
}


quickView <- function(d){
  # Genes
  cat("Top genes\n")
  sel_cv = d$sel_cv
  print(head(sel_cv, 10), digits = 3)

  cat("\nGlobal dNdS\n")
  ## Global dN/dS estimates
  print(d$globaldnds)
  ## annotated table of coding mutations
  # head(d$annotmuts)

  signif_genes = sel_cv[sel_cv$qallsubs_cv<0.1, c("gene_name","qallsubs_cv")]
  rownames(signif_genes) = NULL
  if (nrow(signif_genes)) {
    cat("\nSignificant genes\n")
    print(signif_genes)
  } else cat("\nNo significant genes\n")

}

combineData <- function(snvs, indels){
  snvs <- snvs %>%
    dplyr::select(sample, chrom, pos, ref, alt)
  indels <- indels %>%
    dplyr::select(sample, chrom, pos, ref, alt)
  return(rbind(snvs, indels))
}

add_info <- function(data, attach_info='data/samples_names_conversion.txt', paper=TRUE){
  name_conversion <- read.delim(attach_info, header = F)
  if(paper){
    cat("annotating with publication names\n")
    colnames(name_conversion) <- c("sample_old", "sample_short", "sample", "sex", "assay")
  } else{
    colnames(name_conversion) <- c("sample", "sample_short", "sample_paper", "sex", "assay")
  }

  ann_data <- plyr::join(data, name_conversion, 'sample', type = "left")
  ann_data <- ann_data %>%
    dplyr::mutate(sex = as.factor(ifelse(assay == "whole-gut", "whole-gut", as.character(sex))))
  return(ann_data)
}
