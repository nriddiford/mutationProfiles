#' mutSigs
#'
#' Calculate and plot the mutational signatures accross samples using the package `deconstructSigs`
#' @param samples Calculates and plots mutational signatures on a per-sample basis [Default no]
#' @param pie Plot a pie chart shwoing contribution of each signature to overall profile [Default no] 
#' @import deconstructSigs
#' @import BSgenome
#' @import BSgenome.Dmelanogaster.UCSC.dm6
#' @keywords signatures
#' @export


mutSigs <- function(samples=NA, pie=NA){
  suppressMessages(require(BSgenome.Dmelanogaster.UCSC.dm6))
  suppressMessages(require(deconstructSigs))
  
  if(!exists('scaling_factor')){
    cat("calculationg trinucleotide frequencies in genome\n")
    scaling_factor <-triFreq()
  }
  
  data<-getData()
  genome <- BSgenome.Dmelanogaster.UCSC.dm6
  
  if(is.na(samples)){
    data$tissue = 'All'
    sigs.input <- mut.to.sigs.input(mut.ref = data, sample.id = "tissue", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)
    sig_plot<-whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic, sample.id = 'All',
                              contexts.needed = TRUE,
                              tri.counts.method = scaling_factor
                              )

    cat("Writing to file 'plots/all_signatures.pdf'\n")
    pdf('plots/all_signatures.pdf', width = 20, height = 10)
    plotSignatures(sig_plot)
    dev.off()
    plotSignatures(sig_plot)
    

    if(!is.na(pie)){
      makePie(sig_plot)
    }
  }
  
  else{
  	sigs.input <- mut.to.sigs.input(mut.ref = data, sample.id = "sample", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)
  	cat("sample", "snv_count", sep="\t", "\n")
      for(s in levels(data$sample)) {
        snv_count<-nrow(filter(data, sample == s))
        
        if(snv_count > 50){
          cat(s, snv_count, sep="\t", "\n")
        
          sig_plot<-whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.nature2013, sample.id = s,
                        contexts.needed = TRUE,
                        tri.counts.method = scaling_factor)
          
          outfile<-(paste('plots/', s, '_signatures.pdf', sep = ''))
          cat("Writing to file", outfile, "\n")
          pdf(outfile, width = 20, height = 10)
          plotSignatures(sig_plot)
          dev.off()
          plotSignatures(sig_plot)
          
          if(!is.na(pie)){
            makePie(sig_plot)
          }
        }
      }
  }
}
