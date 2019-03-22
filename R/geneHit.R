#' Print barplot showing samples with SVs affecting specified gene
#' @param infile File to process [Required]
#' @param filter_gene The gene of interest [Default: N]
#' @param plot Show barplot [Default: TRUE]
#' @import stringr ggsci
#' @export
geneHit <- function(..., snv_data = NULL, filter_gene = "N", plot = TRUE) {
  if(missing(snv_data)){
    snv_data<-getData(...)
  }

  geneIn <- function(gene, gene_list) {
    sapply(as.character(gene_list), function(x) gene %in% strsplit(x, ", ")[[1]], USE.NAMES=FALSE)
  }

  gene_hits <- snv_data %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(...,
                  gene == filter_gene)
    dplyr::mutate(n_hits = n()) %>%
    dplyr::select(sample, n_hits, status) %>%
    droplevels()

  sample_names <- snv_data %>%
    # dplyr::filter(...) %>%
    dplyr::group_by(sample) %>%
    dplyr::distinct(sample) %>%
    droplevels()

  sample_names <- levels(sample_names$sample)
  missing_samples = list()

  for(i in 1:length(sample_names)){
    if(!(sample_names[i] %in% levels(gene_hits$sample))){
      missing_samples[i] <- sample_names[i]
    }
  }

  missing_samples <- plyr::compact(missing_samples)

  dat <- data.frame(sample = unlist(missing_samples), n_hits = 0)

  all_samples <- plyr::join(gene_hits, dat, type='full')

  if(plot){
    all_samples$star <- ifelse(all_samples$n_hits==0, "X", '')
    p <- ggplot(all_samples, aes(fct_reorder(sample, -n_hits), n_hits, label = star))
    p <- p + geom_bar(alpha = 0.7, stat = "identity")
    # p <- p + scale_y_continuous("Length (Kb)", expand = c(0, 0.5))
    p <- p + cleanTheme() +
      theme(
        panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
        axis.title.x = element_blank()
      )
    # p <- p + scale_alpha_continuous(range = c(0.1, 1))
    p <- p + ggtitle(paste("Samples with SNVs affecting", filter_gene))
    p <- p + scale_fill_jco()
    p <- p + geom_text(vjust=-5)
    # p <- p + geom_text(data = dat, aes(x=sample_mod, y=50),  label = "X")

    # p <- p + coord_flip()
    # p <- p + scale_y_reverse()

    p
  } else return(gene_hits)
}
