% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneEnrichment.R
\name{geneEnrichment}
\alias{geneEnrichment}
\title{geneEnrichment}
\usage{
geneEnrichment(..., snv_data = NULL,
  gene_lengths_in = "data/gene_lengths.txt", n = 10,
  genome_length = 118274340, write = FALSE)
}
\arguments{
\item{n}{The number of times we need to have seen a gene in our snv_data to view its enrichment score [Default 3]}

\item{genome_length}{The total legnth of the genome [Default 137547960 (chroms 2, 3, 4, X & Y for Drosophila melanogastor Dmel6.12)]}

\item{gene_lengths}{File containing all genes and their lengths (as generated by 'script/genomefeatures.pl') [Default 'data/gene_lengths.txt']}
}
\value{
A snv_data frame with FC scores for all genes seen at least n times in snv snv_data
}
\description{
Calculate the enrichment of snv hits in length-corrected genes
A 'gene_lengths' file must be provided with the following fields (cols 1..6 required)
gene length chrom    start      end      tss scaling_factor
This can be generated using the script 'script/genomic_features.pl' and a genome .gtf file
The defualt genome length is set to the mappable regions of the Drosophila melanogastor Dmel6.12 genome (GEM mappability score > .5)
(118274340). The full, assembled genome legnth for chroms 2/3/4/X/Y is 137547960
}
\details{
Function to calculate fold change enrichment in a set of snv calls correcting for gene length
}
\keyword{enrichment}
