#' mutSigs
#'
#' Calculate and plot the mutational signatures accross samples using the package `deconstructSigs`
#' @param samples Calculates and plots mutational signatures on a per-sample basis [Default no]
#' @param pie Plot a pie chart shwoing contribution of each signature to overall profile [Default no]
#' @import deconstructSigs
#' @import BSgenome.Dmelanogaster.UCSC.dm6
#' @keywords signatures
#' @export

mutSigs <- function(..., snv_data=NULL, by_sample=FALSE, pie=FALSE, write=FALSE){

  if(!exists('scaling_factor')){
    cat("calculating trinucleotide frequencies in genome\n")
    scaling_factor <-triFreq()
  }
  if(missing(snv_data)){
    snv_data<-getData(...)
  } else{
    snv_data <- snv_data %>%
      dplyr::filter(...)
  }

  genome <- BSgenome.Dmelanogaster.UCSC.dm6

  if(!by_sample){
    cat("Plotting for all samples\n")
    snv_data$tissue = 'All'
    sigs.input <- mut.to.sigs.input(mut.ref = snv_data, sample.id = "tissue", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)
    sig_plot<-whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic, sample.id = 'All',
                              contexts.needed = TRUE,
                              tri.counts.method = scaling_factor
    )

    if(write){
      cat("Writing to file 'plots/all_signatures.pdf'\n")
      pdf('plots/all_signatures.pdf', width = 20, height = 10)
      plotSignatures(sig_plot)
      dev.off()
    }
    plotSignatures(sig_plot)

    if(pie){
      makePie(sig_plot)
    }
  }

  else{
    sigs.input <- mut.to.sigs.input(mut.ref = snv_data, sample.id = "sample", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)
    cat("sample", "snv_count", sep="\t", "\n")
    for(s in levels(snv_data$sample)) {
      snv_count<-nrow(filter(snv_data, sample == s))

      if(snv_count > 50){
        cat(s, snv_count, sep="\t", "\n")

        sig_plot<-whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic, sample.id = s,
                                  contexts.needed = TRUE,
                                  tri.counts.method = scaling_factor)

        if(write){
          outfile<-(paste('plots/', s, '_signatures.pdf', sep = ''))
          cat("Writing to file", outfile, "\n")
          pdf(outfile, width = 20, height = 10)
          plotSignatures(sig_plot)
          dev.off()
        }
        plotSignatures(sig_plot)

        if(pie){
          makePie(sig_plot)
        }
      }
    }
  }
}
