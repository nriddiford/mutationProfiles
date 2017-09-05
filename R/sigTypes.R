#' sigTypes
#'
#' Calculate and plot the contribution of mutational signatures accross samples using the package `deconstructSigs`
#' @import deconstructSigs
#' @import data.table
#' @import reshape
#' @import BSgenome.Dmelanogaster.UCSC.dm6
#' @keywords signatures
#' @export

sigTypes <- function(){
  suppressMessages(require(BSgenome.Dmelanogaster.UCSC.dm6))
  suppressMessages(require(deconstructSigs))
  
  if(!exists('scaling_factor')){
    cat("Calculating trinucleotide frequencies in genome\n")
    scaling_factor <- triFreq()
  }
  
  data<-getData()
  genome <- BSgenome.Dmelanogaster.UCSC.dm6
  
  sigs.input <- mut.to.sigs.input(mut.ref = data, sample.id = "sample", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)

  l = list()
  for(s in levels(data$sample)) {
    snv_count<-nrow(filter(data, sample == s))
    
    if(snv_count > 50){

      sig_plot<-whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic, sample.id = s,
                                contexts.needed = TRUE,
                                tri.counts.method = scaling_factor)
      l[[s]] <- sig_plot
    }
  }
  
  mutSigs<-do.call(rbind, l)
  mutSigs<-as.data.frame(mutSigs)

  mutWeights<-mutSigs$weights
  
  mutData<-melt(rbindlist(mutWeights, idcol = 'sample'),
                id = 'sample', variable.name = 'signature', value.name = 'score')
  
  mutData<-filter(mutData, score > 0.1)
  mutData<-droplevels(mutData)
  
  p <- ggplot(mutData[order(mutData$signature),])
  p <- p + geom_bar(aes(reorder(sample, -score), score, fill=signature),colour="black", stat = "identity")
  p <- p + scale_x_discrete("Sample")
  p <- p + scale_y_continuous("Signature contribution", expand = c(0.01, 0.01), breaks=seq(0, 1, by=0.1))
  p <- p + cleanTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  p
}