#' sigTypes
#'
#' Calculate and plot the mutational signatures accross samples using the package `deconstructSigs`
#' @param samples Calculates and plots mutational signatures on a per-sample basis [Default no]
#' @param pie Plot a pie chart shwoing contribution of each signature to overall profile [Default no]
#' @import deconstructSigs
#' @import data.table
#' @import reshape2
#' @import forcats
#' @import BSgenome.Dmelanogaster.UCSC.dm6
#' @keywords signatures
#' @export

sigTypes <- function(..., snv_data=NULL, write=FALSE, signatures=signatures.genome.cosmic.v3.may2019, min_contribution=0.1){

  if(missing(snv_data)){
    snv_data<-getData(...)
  }

  snv_data <- snv_data %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(max(row_number()) > 50) %>%
    droplevels() %>%
    dplyr::ungroup() %>%
    as.data.frame()


  suppressMessages(require(BSgenome.Dmelanogaster.UCSC.dm6))
  suppressMessages(require(deconstructSigs))

  if(!exists('scaling_factor')){
    cat("Calculating trinucleotide frequencies in genome\n")
    scaling_factor <-triFreq()
  }

  genome <- BSgenome.Dmelanogaster.UCSC.dm6

  sigs.input <- mut.to.sigs.input(mut.ref = snv_data, sample.id = "sample", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)

  l = list()
  for(s in levels(snv_data$sample)) {
    snv_count<-nrow(filter(snv_data, sample == s))
    cat(s, snv_count, sep="\t", "\n")

    sig_plot<-whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures, sample.id = s,
                              contexts.needed = TRUE,
                              tri.counts.method = scaling_factor)
    l[[s]] <- sig_plot
  }

  mutSigs<-do.call(rbind, l)
  mutSigs<-as.data.frame(mutSigs)

  mutWeights<-mutSigs$weights

  mutData<-melt(rbindlist(mutWeights, idcol = 'sample'),
                id = 'sample', variable.name = 'signature', value.name = 'score')

  mutData <- mutData %>%
    dplyr::filter(score >= min_contribution) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(total = sum(score)) %>%
    droplevels()

  p <- ggplot(mutData)
  p <- p + geom_bar(aes(fct_reorder(sample, -total), score, fill=signature),colour="black", stat = "identity")
  p <- p + scale_x_discrete("Sample")
  p <- p + scale_y_continuous("Signature contribution", expand = c(0.01, 0.01), breaks=seq(0, 1, by=0.1))
  p <- p + cleanTheme() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
          axis.text = element_text(size=30)
    )

  print(p)

  if(write){
    sigTypes<-paste("sigTypes.pdf")
    cat("Writing file", sigTypes, "\n")
    ggsave(paste("plots/", sigTypes, sep=""), width = 20, height = 10)
  }
  return(mutData)
}


####
# sigTypesPie
####

sigPie <- function() {
  df <- data.frame(
    group = c("Sig3", "Sig5", "Sig8", "Unknown"),
    value = c(21, 14, 25, 40),
    cols = c('#DB8E00', '#64B200', '#00BD5C', '#00BADE'))


  all <- data.frame(
    group = c("Sig3", "Sig8", "Sig9", "Sig21", "Sig25", "Unknown"),
    value = c(29, 17, 10, 7, 7, 30),
    cols = c('#E68613', '#0CB702', '#00BE67', '#ED68ED', '#FF61CC', 'grey'))

  bp <- ggplot(all, aes(x="", y=value, fill = cols)) +
    geom_bar(width = 1, stat = "identity", colour = "white") +
    scale_fill_manual(values = levels(all$cols), labels = levels(all$group))

  pie <- bp + coord_polar("y", start=0)
  pie + cleanTheme() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

}

