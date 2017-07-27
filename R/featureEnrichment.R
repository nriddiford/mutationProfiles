#' featureEnrichment
#'
#' Function to calculate enrichment of snv hits in genomic features
#' @param features File containing total genomic lengths of features [Default 'data/genomic_features.txt']
#' @keywords enrichment
#' @export


featureEnrichment <- function(features='data/genomic_features.txt'){
  genome_features<-read.delim(features, header = T)
  data<-getData()
  breakpoint_count<-nrow(data)
  
  # To condense exon counts into "exon"
  data$feature<-as.factor(gsub("_.*", "", data$feature))
  
  classes_count<-table(data$feature)
  
  class_lengths<-setNames(as.list(genome_features$length), genome_features$feature)
  cat("feature", "observed", "expected", "test", "FC", "sig", "p", "\n")
  
  for (f in levels(data$feature)) {
    feature_fraction<-class_lengths[[f]]/137547960
    feature_expect<-breakpoint_count*(feature_fraction)
    
    if(!is.null(class_lengths[[f]])){
      if(classes_count[f] >= feature_expect){
        stat<-binom.test(x = classes_count[f], n = breakpoint_count, p = feature_fraction, alternative = "greater")
        test<-"enrichment"
      }
      else{
        stat<-binom.test(x = classes_count[f], n = breakpoint_count, p = feature_fraction, alternative = "less")
        test<-"depletion"
      }
      
      ifelse(stat$p.value <= 0.05, sig<-'T', sig<-'F')
      
      p_val<-format.pval(stat$p.value, digits = 3, eps=0.0001)
      fc<-classes_count[f]/feature_expect
      cat(f, classes_count[f], feature_expect, test, fc, sig, p_val, "\n")
      
    }
  } 
}