% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snvinGene.R
\name{snvinGene}
\alias{snvinGene}
\title{snvinGene}
\usage{
snvinGene(..., snv_data = NULL, gene_lengths = "data/gene_lengths.txt",
  gene2plot = "kuz", annotated = TRUE, col_by_status = TRUE,
  write = FALSE)
}
\arguments{
\item{gene_lengths}{File containing all genes and their lengths (as generated by 'script/genomefeatures.pl') [Default 'data/gene_lengths.txt']}

\item{gene2plot}{Name of the gene to plot}
}
\description{
Plot all snvs found in a given gene.
A 'gene_lengths' file must be provided with the following fields (cols 1..6 required)
gene length chrom    start      end      tss scaling_factor
This can be generated using the script 'script/genomic_features.pl' and a genome .gtf file
}
\details{
Plot all snvs found in a given gene
}
\keyword{gene}
