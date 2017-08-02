library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(deconstructSigs)


#' getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords get 
#' @import dplyr
#' @export
#' @return Dataframe

getData <- function(infile = "data/annotated_snvs.txt"){
  data<-read.delim(infile, header = F)
  
  colnames(data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "type", "feature", "gene")
  levels(data$type) <- tolower(levels(data$type))
  #data <- filter(data, type == 'germline')
  
  #filter out samples
  data<-filter(data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  data<-droplevels(data)
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(data)
}


#' cleanTheme
#'
#' Clean theme for plotting
#' @param base_size Base font size [Default 12]
#' @import ggplot2
#' @keywords theme
#' @export

cleanTheme <- function(base_size = 12){
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text = element_text(size=12),
    axis.title = element_text(size=30)
  )
}


#' genTris
#'
#' This function returns all possible trinucleotide combinations
#' @keywords trinucleotides
#' @export
#' @return Character string containing all 96 trinucleotides
#' genTris()

genTris <- function(){
  all.tri = c()
    for(i in c("A", "C", "G", "T")){
      for(j in c("C", "T")){
        for(k in c("A", "C", "G", "T")){
          if(j != k){
            for(l in c("A", "C", "G", "T")){
              tmp = paste(i, "[", j, ">", k, "]", l, sep = "")
              all.tri = c(all.tri, tmp)
            }
          }
        }
      }
    }
  all.tri <- all.tri[order(substr(all.tri, 3, 5))]
  return(all.tri)
}


#' setCols
#'
#' Show top hit genes
#' @import RColorBrewer
#' @param df Dataframe [Required]
#' @param col Column of dataframe. Colours will be set to levels(df$cols) [Required]
#' @keywords cols
#' @export

setCols <- function(df, col){
  names<-levels(df[[col]])
  cat("Setting colour levles:", names, "\n")
  level_number<-length(names)
  mycols<-brewer.pal(level_number, "Set2")
  names(mycols) <- names
  colScale <- scale_fill_manual(name = col,values = mycols)
  return(colScale)
}


#' chromDist
#'
#' Plot genome-wide snv distribution 
#' @import ggplot2
#' @keywords distribution
#' @export

chromDist <- function(object=NA, notch=0){
  data<-getData()
  ext<-'.pdf'
  if(is.na(object)){
    object<-'grouped_trans'
    cols<-setCols(data, "grouped_trans")
  }
  
  if(notch){
    data<-exclude_notch()
    ext<-'_excl.N.pdf'
  }

  cat("Plotting snvs by", object, "\n")
  
  p<-ggplot(data)
  p<-p + geom_histogram(aes(pos/1000000, fill = get(object)), binwidth=0.1, alpha = 0.8)  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33),expand = c(0.01, 0.01))
  p<-p + scale_y_continuous("Number of snvs", expand = c(0.01, 0.01))
  p<-p + cleanTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15)
    )
  
  if (object == 'grouped_trans'){
    p<-p + cols
  }  
  chrom_outfile<-paste("snv_dist_genome_by_", object, ext, sep = "")
  cat("Writing file", chrom_outfile, "\n")
  ggsave(paste("plots/", chrom_outfile, sep=""), width = 20, height = 10)
  
  p
}


#' featuresHit
#'
#' Show top hit features
#' @import ggplot2
#' @keywords features
#' @export

featuresHit <- function(){
  data<-getData()
  
  # To condense exon counts into "exon"
  data$feature<-as.factor(gsub("exon_.*", "exon", data$feature))
  
  # Reoders descending
  data$feature<-factor(data$feature, levels = names(sort(table(data$feature), decreasing = TRUE)))
  
  #cols<-set_cols(data, "feature")
  
  p<-ggplot(data)
  p<-p + geom_bar(aes(feature, fill = feature))
  #p<-p + cols
  p<-p + cleanTheme() +
    theme(axis.title.x=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous(expand = c(0.01, 0.01))
  
  features_outfile<-paste("hit_features_count.pdf")
  cat("Writing file", features_outfile, "\n")
  ggsave(paste("plots/", features_outfile, sep=""), width = 20, height = 10)
  
  p
}


#' geneHit
#'
#' Show top hit genes
#' @import dplyr
#' @keywords gene
#' @param n Show top n hits [Default 10] 
#' @export

geneHit <- function(n=10){
  data<-getData()
  data<-filter(data, gene != "intergenic")
  
  hit_count<-as.data.frame(sort(table(unlist(data$gene)), decreasing = T))
  
  colnames(hit_count)<- c("gene", "count")
  head(hit_count, n)
}


#' genomeSnvs
#'
#' Plot snvs accross genome
#' @import ggplot2
#' @keywords genome
#' @export

genomeSnvs <- function(){
  data<-getData()
  data<-filter(data, chrom != "4")
  p<-ggplot(data)
  p<-p + geom_point(aes(pos/1000000, sample, colour = trans))
  # p<-p + guides(color = FALSE)
  p<-p + theme(axis.text.x = element_text(angle=45, hjust = 1))
  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  
  p
}


#' mutSigs
#'
#' Calculate and plot the mutational signatures accross samples using the package `deconstructSigs`
#' @param samples Calculates and plots mutational signatures on a per-sample basis [Default no]
#' @param pie Plot a pie chart shwoing contribution of each signature to overall profile [Default no] 
#' @import deconstructSigs
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


#' mutSpectrum
#'
#' Plots the mutations spectrum for all samples combined
#' @import ggplot2
#' @keywords spectrum
#' @export

mutSpectrum <- function(){
  data<-getData()
  cat("Showing global contribution of tri class to mutation load", "\n")
 
  p<-ggplot(data)
  p<-p + geom_bar(aes(x = decomposed_tri, y = (..count..)/sum(..count..), group = decomposed_tri, fill = grouped_trans), position="dodge",stat="count")
  p<-p + scale_y_continuous("Relative contribution to mutation load", expand = c(0.0, .0005))
  p<-p + scale_x_discrete("Genomic context", expand = c(.005, .005))
  p<-p + cleanTheme() + 
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15)
          )
  p<-p + labs(fill="Mutation class")
  p<-p + facet_wrap(~grouped_trans, ncol = 3, scale = "free_x" )
  
  mut_spectrum<-paste("mutation_spectrum.pdf")
  cat("Writing file", mut_spectrum, "\n")
  ggsave(paste("plots/", mut_spectrum, sep=""), width = 20, height = 10)
  p
}


#' snvinGene
#'
#' Plot all snvs found in a given gene
#' @param gene_lengths File containing all genes and their lengths (as generated by 'script/genomefeatures.pl') [Default 'data/gene_lengths.txt']
#' @param gene2plot Name of the gene to plot
#' @import ggplot2
#' @keywords gene
#' @export

snvinGene <- function(gene_lengths="data/gene_lengths.txt", gene2plot='dnc'){
  gene_lengths<-read.delim(gene_lengths, header = T)
  region<-filter(gene_lengths, gene == gene2plot)
  
  wStart<-(region$start - 10000)
  wEnd<-(region$end + 10000)
  wChrom<-region$chrom
  wTss<-suppressWarnings(as.numeric(levels(region$tss))[region$tss])
  data<-getData()
  data<-filter(data, chrom == wChrom & pos >= wStart & pos <= wEnd)
  
  if(nrow(data) == 0){
    stop(paste("There are no snvs in", gene2plot, "- Exiting", "\n"))
  }
  
  p<-ggplot(data)
  p<-p + geom_point(aes(pos/1000000, sample, colour = trans, size = 1.5), position=position_jitter(width=0.00005, height=0))
  p<-p + guides(size = FALSE, sample = FALSE)
  p<-p + cleanTheme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(round(wStart/1000000, digits = 2),round(wEnd/1000000, digits = 2),by=0.05), limits=c(wStart/1000000, wEnd/1000000))
  p<-p + annotate("rect", xmin=region$start/1000000, xmax=region$end/1000000, ymin=0, ymax=0.1, alpha=.2, fill="skyblue")

  p<-p + geom_segment(aes(x = wTss/1000000, y = 0, xend= wTss/1000000, yend = 0.1), colour="red")
  middle<-((wEnd/1000000+wStart/1000000)/2)
  p <- p + annotate("text", x = middle, y = 0.05, label=gene2plot, size=6)
  p<-p + ggtitle(paste("Chromosome:", wChrom))
  
  p
}


#' tssDist
#'
#' Plot distance to TSS distribution
#' @param gene_lengths File containing all genes and their lengths (as generated by 'script/genomefeatures.pl') [Default 'data/gene_lengths.txt']
#' @import ggplot2
#' @keywords tss
#' @export

tssDist <- function(gene_lengths="data/gene_lengths.txt"){
  gene_lengths<-read.delim(gene_lengths, header = T)
  
  
  data$dist2tss<-
  
  
  data<-getData()
  data<-filter(data, chrom == wChrom & pos >= wStart & pos <= wEnd)
  
  if(nrow(data) == 0){
    stop(paste("There are no snvs in", gene2plot, "- Exiting", "\n"))
  }
  
  p<-ggplot(data)
  p<-p + geom_point(aes(pos/1000000, sample, colour = trans, size = 1.5), position=position_jitter(width=0.00005, height=0))
  p<-p + guides(size = FALSE, sample = FALSE)
  p<-p + cleanTheme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(round(wStart/1000000, digits = 2),round(wEnd/1000000, digits = 2),by=0.05), limits=c(wStart/1000000, wEnd/1000000))
  p<-p + annotate("rect", xmin=region$start/1000000, xmax=region$end/1000000, ymin=0, ymax=0.1, alpha=.2, fill="skyblue")
  
  middle<-((wEnd/1000000+wStart/1000000)/2)
  p <- p + annotate("text", x = middle, y = 0.05, label=gene2plot, size=6)
  p<-p + ggtitle(paste("Chromosome:", wChrom))
  
  p
}




#' samplesPlot
#'
#' Plot the snv distribution for each sample
#' @import ggplot2
#' @param count Output total counts instead of frequency if set [Default no] 
#' @keywords spectrum
#' @export

samplesPlot <- function(count=NA){
  data<-getData()

  mut_class<-c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

  p<-ggplot(data)
  
  if(is.na(count)){
    p<-p + geom_bar(aes(x = grouped_trans, y = (..count..)/sum(..count..), group = sample, fill = sample), position="dodge",stat="count")
    p<-p + scale_y_continuous("Relative contribution to total mutation load", expand = c(0.0, .001))
    tag='_freq'
  }
  else{
    p<-p + geom_bar(aes(x = grouped_trans, y = ..count.., group = sample, fill = sample), position="dodge",stat="count")
    p<-p + scale_y_continuous("Count", expand = c(0.0, .001))
    tag='_count'
  }
  p<-p + scale_x_discrete("Mutation class", limits=mut_class)
  p<-p + cleanTheme() + 
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 10)
          )
  p<-p + facet_wrap(~sample, ncol = 4, scale = "free_x" )
  
  samples_mut_spect<-paste("mutation_spectrum_samples", tag, ".pdf", sep = '')
  cat("Writing file", samples_mut_spect, "\n")
  ggsave(paste("plots/", samples_mut_spect, sep=""), width = 20, height = 10)
  p
}


#' snvStats
#'
#' Calculate some basic stats for snv data
#' @import dplyr
#' @keywords stats
#' @export

snvStats <- function(){
  data<-getData()
  cat("sample", "snvs", sep='\t', "\n")
  rank<-sort(table(data$sample), decreasing = TRUE)
  rank<-as.array(rank)
  
  for (i in 1:nrow(rank)){
    cat(names(rank[i]), rank[i], sep='\t', "\n")
  }
  
  all_ts<-nrow(filter(data, trans == "A>G" | trans == "C>T" | trans == "G>A" | trans == "T>C"))
  all_tv<-nrow(filter(data, trans != "A>G" & trans != "C>T" & trans != "G>A" & trans != "T>C"))
  ts_tv<-all_ts/all_tv
  cat("ts/tv =", ts_tv)
}


#' triFreq
#'
#' This function counts the number of times each triunucleotide is found in a supplied genome
#' @param genome BS.genome file defaults to BSgenome.Dmelanogaster.UCSC.dm6
#' @param count Output total counts instead of frequency if set [Default no] 
#' @import dplyr
#' @keywords trinucleotides
#' @export
#' @return Dataframe of trinucs and freqs (or counts if count=1)

triFreq <- function(genome=NA, count=NA){
  if(is.na(genome)){
    cat("No genome specfied, defaulting to 'BSgenome.Dmelanogaster.UCSC.dm6'\n")
    library(BSgenome.Dmelanogaster.UCSC.dm6, quietly = TRUE)
    genome <- BSgenome.Dmelanogaster.UCSC.dm6
  }
  
  params <- new("BSParams", X = Dmelanogaster, FUN = trinucleotideFrequency, exclude = c("M", "_"), simplify = TRUE)
  data<-as.data.frame(bsapply(params))
  data$genome<-as.integer(rowSums(data))
  data$genome_adj<-(data$genome*2)
  
  if(!is.na(count)){
    tri_count<-data['genome_adj']
    tri_count<-cbind(tri = rownames(tri_count), tri_count)
    colnames(tri_count) <- c("tri", "count")
    rownames(tri_count) <- NULL
    return(tri_count)
  }
  else{
    data$x <- (1/data$genome)
    scaling_factor<-data['x']
    return(scaling_factor)
  }
  
}


#' featureEnrichment
#'
#' Function to calculate enrichment of snv hits in genomic features
#' @param features File containing total genomic lengths of features [Default 'data/genomic_features.txt']
#' @param genome_length The total legnth of the genome [Default 137547960 (chroms 2, 3, 4, X & Y for Drosophila melanogastor Dmel6.12)]
#' @keywords enrichment
#' @import dplyr
#' @return A data frame with FC scores for all genes seen at least n times in snv data
#' @export 

featureEnrichment <- function(features='data/genomic_features.txt', genome_length=137547960){
  genome_features<-read.delim(features, header = T)
  data<-getData()
  mutCount<-nrow(data)
  
  # To condense exon counts into "exon"
  data$feature<-as.factor(gsub("exon_.*", "exon", data$feature))
  
  classCount<-table(data$feature)
  classLengths<-setNames(as.list(genome_features$length), genome_features$feature)

  fun <- function(f) {
    # Calculate the fraction of geneome occupied by each feature
    featureFraction<-classLengths[[f]]/genome_length
    
    # How many times should we expect to see this feature hit in our data (given number of obs. and fraction)?
    featureExpect<-(mutCount*featureFraction)
    
    # observed/expected 
    fc<-classCount[[f]]/featureExpect
    fc<-round(fc,digits=1)
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
      sig_val<-'F'
      if(stat$p.value <= 0.05){ sig_val<-'T'}
      p_val<-format.pval(stat$p.value, digits = 3, eps=0.0001)
      list(feature = f, observed = classCount[f], expected = featureExpect, fc = fc, test = test, sig = sig_val, p_val = p_val)
    }
  }
  
  enriched<-lapply(levels(data$feature), fun)
  enriched<-do.call(rbind, enriched)
  featuresFC<-as.data.frame(enriched)
  # Sort by FC value
  featuresFC<-arrange(featuresFC,desc(as.integer(fc)))
  return(featuresFC)
}


#' geneEnrichment
#'
#' Function to calculate fold change enrichment in a set of snv calls correcting for gene length
#' @param gene_lengths File containing all genes and their lengths (as generated by 'script/genomefeatures.pl') [Default 'data/gene_lengths.txt']
#' @param n The number of times we need to have seen a gene in our data to view its enrichment score [Default 2]
#' @param genome_length The total legnth of the genome [Default 137547960 (chroms 2, 3, 4, X & Y for Drosophila melanogastor Dmel6.12)]
#' @keywords enrichment
#' @import dplyr
#' @return A data frame with FC scores for all genes seen at least n times in snv data
#' @export 

geneEnrichment <- function(gene_lengths="data/gene_lengths.txt", n=2, genome_length=137547960){
  gene_lengths<-read.delim(gene_lengths, header = T)
  data<-getData()
  data<-filter(data, gene != "intergenic")
  
  snv_count<-nrow(data)
  
  hit_genes<-table(data$gene)
  genes<-setNames(as.list(gene_lengths$length), gene_lengths$gene)
  
  fun <- function(g) {
    # Calculate the fraction of geneome occupied by each gene
    genefraction<-genes[[g]]/genome_length
    
    # How many times should we expect to see this gene hit in our data (given number of obs. and fraction)?
    gene_expect<-snv_count*(genefraction)
    
    # observed/expected 
    fc<-hit_genes[[g]]/gene_expect
    fc<-round(fc,digits=1)
    gene_expect<-round(gene_expect,digits=3)
    list(gene = g, length = genes[[g]], observed = hit_genes[g], expected = gene_expect, fc = fc)
  }
  
  enriched<-lapply(levels(data$gene), fun)
  enriched<-do.call(rbind, enriched)
  genesFC<-as.data.frame(enriched)
  # Filter for genes with few observations
  genesFC<-filter(genesFC, observed > n)
  # Sort by FC value
  genesFC<-arrange(genesFC,desc(as.integer(fc)))
  return(genesFC)
}




