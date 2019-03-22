
genes_of_interest <- function(df, annotation_file = '~/Desktop/gene2bed/merged_groups.bed', plot=F){
  annotation <- read.delim(annotation_file, header = F)
  colnames(annotation) <- c("chrom", "start", "end", "gene", "id", "groups")

  annotation <- annotation %>%
    dplyr::select(gene, groups)

  the_specials <- dplyr::inner_join(df, annotation, by="gene")

  if(plot){
  p <- ggplot(the_specials)
  p <- p + geom_bar(aes(sample, gene))
  }

}
