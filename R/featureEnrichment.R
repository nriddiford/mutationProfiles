#' featureEnrichment
#'
#' Function to calculate enrichment of snv hits in genomic features
#' @description Calculate the enrichment of snv hits in genomic features
#' A 'features' file must be provided with the follwing format:
#' feature	length	percentage
#' This can be generated using the script 'script/genomic_features.pl' and a genome .gtf file
#' The defualt genome length is set to the mappable regions of the Drosophila melanogastor Dmel6.12 genome (GEM mappability score > .5)
#' (118274340). The full, assembled genome legnth for chroms 2/3/4/X/Y is 137547960
#' @param features File containing total genomic lengths of features [Default 'data/genomic_features.txt']
#' @param genome_length The total legnth of the genome [Default 118274340 (mappable regions on chroms 2, 3, 4, X & Y for Drosophila melanogastor Dmel6.12)]
#' @keywords enrichment
#' @import dplyr ggpubr
#' @return A snv_data frame with FC scores for all genes seen at least n times in snv snv_data
#' @export
featureEnrichment <- function(..., features='data/genomic_features.txt', genome_length=118274340, write=FALSE){
  genome_features<-read.delim(features, header = T)
  snv_data<-getData(...)
  mutCount<-nrow(snv_data)

  # To condense exon counts into "exon"
  snv_data$feature<-as.factor(gsub("exon_.*", "exon", snv_data$feature))

  classCount<-table(snv_data$feature)
  classLengths<-setNames(as.list(genome_features$length), genome_features$feature)

  fun <- function(f) {
    # Calculate the fraction of geneome occupied by each feature
    featureFraction<-classLengths[[f]]/genome_length

    # How many times should we expect to see this feature hit in our snv_data (given number of obs. and fraction)?
    featureExpect<-(mutCount*featureFraction)

    # observed/expected
    fc<-classCount[[f]]/featureExpect
    Log2FC<-log2(fc)
    featureExpect<-round(featureExpect,digits=3)

    # Binomial test
    if(!is.null(classLengths[[f]])){
      if(classCount[f] >= featureExpect){
        stat<-binom.test(x = classCount[f], n = mutCount, p = featureFraction, alternative = "greater")
        test<-"enrichment"
      }
      else{
        stat<-binom.test(x = classCount[f], n = mutCount, p = featureFraction, alternative = "less")
        test<-"depletion"
      }
      sig_val <- ifelse(stat$p.value <= 0.001, "***",
                        ifelse(stat$p.value <= 0.01, "**",
                               ifelse(stat$p.value <= 0.05, "*", "")))

      p_val<-format.pval(stat$p.value, digits = 3, eps=0.0001)
      list(feature = f, observed = classCount[f], expected = featureExpect, Log2FC = Log2FC, test = test, sig = sig_val, p_val = p_val)
    }
  }

  enriched<-lapply(levels(snv_data$feature), fun)
  enriched<-do.call(rbind, enriched)
  featuresFC<-as.data.frame(enriched)
  # Sort by FC value
  featuresFC<-dplyr::arrange(featuresFC,desc(abs(as.numeric(Log2FC))))
  featuresFC$Log2FC<-round(as.numeric(featuresFC$Log2FC), 1)

  if(write){
    featuresFC <- filter(featuresFC, observed >= 5)
    first.step <- lapply(featuresFC, unlist)
    second.step <- as.data.frame(first.step, stringsAsFactors = F)
    ggpubr::ggtexttable(second.step, rows = NULL, theme = ttheme("mGreen"))
    feat_enrichment_table <- paste("feature_enrichment_table.tiff")
    cat("Writing to file: ", 'plots/', feat_enrichment_table, sep = '')
    ggsave(paste("plots/", feat_enrichment_table, sep=""), width = 5.5, height = (nrow(featuresFC)/3), dpi=300)
  } else{
    return(featuresFC)
  }
}


featureEnrichmentPlot <- function(write=FALSE) {
  feature_enrichment<-featureEnrichment()

  feature_enrichment$feature <- as.character(feature_enrichment$feature)
  feature_enrichment$Log2FC <- as.numeric(feature_enrichment$Log2FC)

  feature_enrichment <- transform(feature_enrichment, feature = reorder(feature, -Log2FC))

  feature_enrichment <- filter(feature_enrichment, observed >= 5)

  # Custom sorting
  # feature_enrichment$feature <- factor(feature_enrichment$feature, levels=c("intron", "intergenic", "exon", "3UTR", "ncRNA", "5UTR"))

  p<-ggplot(feature_enrichment)
  p<-p + geom_bar(aes(feature, Log2FC, fill = as.character(test)), stat="identity")
  p<-p + guides(fill=FALSE)
  p<-p + ylim(-2,2)
  p<-p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=20)
    )

  if(write){
    feat_plot <- paste("feat_plot.pdf")
    cat("Writing file", feat_plot, "\n")
    ggsave(paste("plots/", feat_plot, sep=""), width = 5, height = 10)
  }
  p
}
