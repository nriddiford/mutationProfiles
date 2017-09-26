list.of.packages <- c('ggplot2', 'dplyr', 'plyr', 'RColorBrewer',
                      'BSgenome.Dmelanogaster.UCSC.dm6', 'deconstructSigs',
                      'reshape', 'data.table', 'ggpubr', 'plotly', 'grid', 'VennDiagram')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat('Installing missing packages...\n')
  install.packages(new.packages)
}
cat('Silently loading packages...')
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))
suppressMessages(library(deconstructSigs))
suppressMessages(library(reshape))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(plotly))
suppressMessages(library(grid))
suppressMessages(library(VennDiagram))

set.seed(42)


#' getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords get 
#' @import dplyr
#' @export
#' @return Dataframe

getData <- function(infile = "data/annotated_snvs.txt", expression_data='data/isc_genes_rnaSeq.csv'){
  snv_data<-read.delim(infile, header = F)
  colnames(snv_data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "a_freq", "caller", "feature", "gene", "id")

  # Read in tissue specific expression data
  seq_data<-read.csv(header = F, expression_data)
  colnames(seq_data)<-c('fpkm', 'id')
  
  # Left might be better(...)
  snv_data <- join(snv_data,seq_data,"id", type = 'left')
  
  snv_data$fpkm<-round(snv_data$fpkm, 1)
  
  # Order by FPKM
  snv_data<- arrange(snv_data, desc(fpkm))
  
  # Filter for genes expressed in RNA-Seq data
  #snv_data<-filter(snv_data, !is.na(fpkm) & fpkm > 0.1)

  # Filter on allele freq
  #snv_data<-filter(snv_data, is.na(a_freq))
 # snv_data<-filter(snv_data, a_freq >= 0.20)
  
  # Filter out samples
  snv_data<-filter(snv_data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  #snv_data<-filter(snv_data, sample != "A373R11" & sample != 'A373R13')
  
  # Find vars called by both Mu and Var 
  # Must also filter one of these calls out...
  snv_data$dups<-duplicated(snv_data[,1:3])
  snv_data<-mutate(snv_data, caller = ifelse(dups == "TRUE", 'varscan2_mutect2' , as.character(caller)))
  
  # Filter for calls made by both V and M
  # snv_data<-filter(snv_data, caller == 'varscan2_mutect2')
  
  # Filter for old/new data
  #snv_data <- filter(snv_data, grepl("^A|H", sample))
  
  snv_data<-droplevels(snv_data)
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(snv_data)
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
#' Get colours for n levels
#' @import RColorBrewer
#' @param df Dataframe [Required]
#' @param col Column of snv_dataframe. Colours will be set to levels(df$cols) [Required]
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


#' snvStats
#'
#' Calculate some basic stats for snv snv_data
#' @import dplyr
#' @keywords stats
#' @export

snvStats <- function(){
  snv_data<-getData()
  cat("sample", "snvs", sep='\t', "\n")
  rank<-sort(table(snv_data$sample), decreasing = TRUE)
  rank<-as.array(rank)
  
  total=0
  
  scores=list()
  for (i in 1:nrow(rank)){
    cat(names(rank[i]), rank[i], sep='\t', "\n")
    total<-total + rank[i]
    scores[i]<-rank[i]
  }
  cat('--------------', '\n')
  scores<-unlist(scores)
  
  mean<-as.integer(mean(scores))
  med<-as.integer(median(scores))
  
  cat('total', total, sep='\t', '\n')
  cat('samples', nrow(rank), sep='\t', '\n')
  
  cat('--------------', '\n')
  cat('mean', mean, sep='\t', '\n')
  cat('median', med, sep='\t', '\n')
  
  cat('\n')
  all_ts<-nrow(filter(snv_data, trans == "A>G" | trans == "C>T" | trans == "G>A" | trans == "T>C"))
  all_tv<-nrow(filter(snv_data, trans != "A>G" & trans != "C>T" & trans != "G>A" & trans != "T>C"))
  ts_tv<-round((all_ts/all_tv), digits=2)
  cat("ts/tv = ", ts_tv,  sep='', '\n')

}



#' rainfall
#'
#' Plot log10 distances between snvs as rainfall plot
#' @import ggplot2
#' @keywords rainfall
#' @export

rainfall <- function(){
  snv_data<-getData()
  
  snv_data<-filter(snv_data, sample != "A373R11" & sample != 'A373R13')
  distances<-do.call(rbind, lapply(split(snv_data[order(snv_data$chrom, snv_data$pos),], snv_data$chrom[order(snv_data$chrom, snv_data$pos)]),
                                   function(a) 
                                     data.frame(a,
                                                dist=c(diff(a$pos), NA),
                                                logdist = c(log10(diff(a$pos)), NA))
  )
  )
  
  
  distances$logdist[is.infinite(distances$logdist)] <- 0
  distances<-filter(distances, chrom != 4)
  
  p<-ggplot(distances)
  p<-p + geom_point(aes(pos/1000000, logdist, colour = grouped_trans))
  p <- p + cleanTheme() +
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          strip.text = element_text(size=20)
    )
  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 6)
  #p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,max(distances$pos),by=10))
  
  rainfall_out<-paste("rainfall.pdf")
  cat("Writing file", rainfall_out, "\n")
  ggsave(paste("plots/", rainfall_out, sep=""), width = 20, height = 5)
  
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
  snv_data<-getData()
  
  mut_class<-c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  
  p<-ggplot(snv_data)
  
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

#calledSnvs

calledSnvs <- function(){
  snv_data<-getData()
  calls<-table(snv_data$caller)
  calls<-as.data.frame(unlist(calls))
  calls$Var1 <- as.factor(calls$Var1)
  
  
  grid.newpage()
  draw.pairwise.venn(area1 = calls$Freq[calls$Var1 == 'mutect2'],
                     area2 = calls$Freq[calls$Var1 == 'varscan2'],
                     cross.area = calls$Freq[calls$Var1 == 'varscan2_mutect2'],
                     category = c("Mutect2","Varscan2"),
                     #lty = rep('blank', 2),
                     lwd = rep(0.3, 2), 
                     cex = rep(2, 3),
                     cat.cex = rep(2, 2),
                     fill = c("#E7B800", "#00AFBB"),
                     alpha = rep(0.4, 2),
                     cat.pos = c(0, 0),
                     #cat.dist = rep(0.025, 2)
                     ext.text = 'FALSE'
  )
  
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
  
  snv_data<-getData()
  genome <- BSgenome.Dmelanogaster.UCSC.dm6
  
  if(is.na(samples)){
    snv_data$tissue = 'All'
    sigs.input <- mut.to.sigs.input(mut.ref = snv_data, sample.id = "tissue", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)
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
    sigs.input <- mut.to.sigs.input(mut.ref = snv_data, sample.id = "sample", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)
    cat("sample", "snv_count", sep="\t", "\n")
    for(s in levels(snv_data$sample)) {
      snv_count<-nrow(filter(snv_data, sample == s))
      
      if(snv_count > 50){
        cat(s, snv_count, sep="\t", "\n")
        
        sig_plot<-whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic, sample.id = s,
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


#' sigTypes
#'
#' Calculate and plot the mutational signatures accross samples using the package `deconstructSigs`
#' @param samples Calculates and plots mutational signatures on a per-sample basis [Default no]
#' @param pie Plot a pie chart shwoing contribution of each signature to overall profile [Default no] 
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
    scaling_factor <-triFreq()
  }
  
  snv_data<-getData()
  genome <- BSgenome.Dmelanogaster.UCSC.dm6
  
  sigs.input <- mut.to.sigs.input(mut.ref = snv_data, sample.id = "sample", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)
  
  l = list()
  for(s in levels(snv_data$sample)) {
    snv_count<-nrow(filter(snv_data, sample == s))
    
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
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=30)
    )
  
  sigTypes<-paste("sigTypes.pdf")
  cat("Writing file", sigTypes, "\n")
  ggsave(paste("plots/", sigTypes, sep=""), width = 20, height = 10)
  
  p
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





#' mutSpectrum
#'
#' Plots the mutations spectrum for all samples combined
#' @import ggplot2
#' @keywords spectrum
#' @export

mutSpectrum <- function(){
  snv_data<-getData()
  cat("Showing global contribution of tri class to mutation load", "\n")
  
  p<-ggplot(snv_data)
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
#' @import dplyr
#' @return A snv_data frame with FC scores for all genes seen at least n times in snv snv_data
#' @export 

featureEnrichment <- function(features='data/genomic_features.txt', genome_length=118274340, print=NA){
  genome_features<-read.delim(features, header = T)
  snv_data<-getData()
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
  
  enriched<-lapply(levels(snv_data$feature), fun)
  enriched<-do.call(rbind, enriched)
  featuresFC<-as.data.frame(enriched)
  # Sort by FC value
  featuresFC<-arrange(featuresFC,desc(as.integer(fc)))
  
  if(!is.na(print)){
    first.step <- lapply(featuresFC, unlist) 
    second.step <- as.data.frame(first.step, stringsAsFactors = F)
    
    ggtexttable(second.step, rows = NULL, theme = ttheme("mBlue"))
    
    feat_enrichment_table <- paste("feature_enrichment_table.pdf")
    cat("Writing to file: ", 'plots/', feat_enrichment_table, sep = '')
    
    ggsave(paste("plots/", feat_enrichment_table, sep=""), width = 5, height = (nrow(featuresFC)/3))
    
  }
  
  else{ return(featuresFC) }
}


featureEnrichmentPlot <- function() {
  feature_enrichment<-featureEnrichment()
  
  feature_enrichment$Log2FC <- log2(as.numeric(feature_enrichment$fc))
  
  feature_enrichment$feature <- as.character(feature_enrichment$feature)
  feature_enrichment$fc <- as.numeric(feature_enrichment$fc)
  
  feature_enrichment <- transform(feature_enrichment, feature = reorder(feature, -fc))
  
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
  
  feat_plot <- paste("feat_plot.pdf")
  cat("Writing file", feat_plot, "\n")
  ggsave(paste("plots/", feat_plot, sep=""), width = 5, height = 10)
  p
  
}


#' geneEnrichment
#'
#' Function to calculate fold change enrichment in a set of snv calls correcting for gene length
#' @description Calculate the enrichment of snv hits in length-corrected genes
#' A 'gene_lengths' file must be provided with the following fields (cols 1..6 required)
#' gene length chrom    start      end      tss scaling_factor
#' This can be generated using the script 'script/genomic_features.pl' and a genome .gtf file
#' The defualt genome length is set to the mappable regions of the Drosophila melanogastor Dmel6.12 genome (GEM mappability score > .5)
#' (118274340). The full, assembled genome legnth for chroms 2/3/4/X/Y is 137547960
#' @param gene_lengths File containing all genes and their lengths (as generated by 'script/genomefeatures.pl') [Default 'data/gene_lengths.txt']
#' @param n The number of times we need to have seen a gene in our snv_data to view its enrichment score [Default 3]
#' @param genome_length The total legnth of the genome [Default 137547960 (chroms 2, 3, 4, X & Y for Drosophila melanogastor Dmel6.12)]
#' @keywords enrichment
#' @import dplyr
#' @import ggpubr
#' @return A snv_data frame with FC scores for all genes seen at least n times in snv snv_data
#' @export 

geneEnrichment <- function(gene_lengths_in="data/gene_lengths.txt", n=7, genome_length=118274340, print=NA){
  gene_lengths<-read.delim(gene_lengths_in, header = T)
  snv_data<-getData()
  snv_data<-filter(snv_data, gene != "intergenic")
  
  gene_lengths<-filter(gene_lengths, length > 1000)
  gene_lengths<-droplevels(gene_lengths)
  
  snv_data<-droplevels(snv_data)
  
  snv_count<-nrow(snv_data)
  
  
  genes<-setNames(as.list(gene_lengths$length), gene_lengths$gene)
  
  # Only keep relevant cols
  gene_lengths<-gene_lengths[,c("gene","length")]
  gene_lengths<-droplevels(gene_lengths)
 
  snv_data<-join(gene_lengths, snv_data, 'gene', type = 'left')
  
  snv_data$fpkm <- ifelse(snv_data$fpkm=='NULL' | snv_data$fpkm=='NA' | is.na(snv_data$fpkm), 0, snv_data$fpkm)
  snv_data$observed <- ifelse(is.numeric(snv_data$observed), snv_data$observed, 0)
 
  hit_genes<-table(factor(snv_data$gene, levels = levels(snv_data$gene) ))
  expression<-setNames(as.list(snv_data$fpkm), snv_data$gene)
  
  fun <- function(g) {
    # Calculate the fraction of geneome occupied by each gene
    genefraction<-genes[[g]]/genome_length
    
    # How many times should we expect to see this gene hit in our snv_data (given number of obs. and fraction of genome)?
    gene_expect<-snv_count*(genefraction)
    
    # observed/expected 
    fc<-hit_genes[[g]]/gene_expect
    log2FC = log2(fc)
    
    gene_expect<-round(gene_expect,digits=3)
    list(gene = g, length = genes[[g]], fpkm = expression[[g]],  observed = hit_genes[g], expected = gene_expect, fc = fc, log2FC = log2FC)
  }
  
  enriched<-lapply(levels(snv_data$gene), fun)
  enriched<-do.call(rbind, enriched)
  genesFC<-as.data.frame(enriched)
  # Filter for genes with few observations
  genesFC<-filter(genesFC, observed >= n)
  genesFC<-droplevels(genesFC)
  genesFC$expected<-round(as.numeric(genesFC$expected),digits=3)
  genesFC$fc<-round(as.numeric(genesFC$fc), 2)
  genesFC$log2FC<-round(as.numeric(genesFC$log2FC), 2)
  # Sort by FC value
  genesFC<-arrange(genesFC,desc(fc))

  
  
  if(!is.na(print)){
    cat("printing")
    first.step <- lapply(genesFC, unlist) 
    second.step <- as.data.frame(first.step, stringsAsFactors = F)
    arrange(second.step,desc(as.integer(log2FC)))
    
    ggtexttable(second.step, rows = NULL, theme = ttheme("mOrange"))
    
    gene_enrichment_table <- paste("gene_enrichment_table.pdf")
    ggsave(paste("plots/", gene_enrichment_table, sep=""), width = 5.2, height = (nrow(genesFC)/3))
  }
  
  else{ return(genesFC) }
}


geneEnrichmentPlot <- function(n=0) {
  gene_enrichment<-geneEnrichment(n=n)
  
  #gene_enrichment<-filter(gene_enrichment, fpkm > 0)
  
  gene_enrichment$gene <- as.character(gene_enrichment$gene)
  gene_enrichment$fc <- as.numeric(gene_enrichment$fc)
  gene_enrichment$log2FC <- as.numeric(gene_enrichment$log2FC)
  
  gene_enrichment <- transform(gene_enrichment, gene = reorder(gene, -fc))
  
  gene_enrichment$test <- ifelse(gene_enrichment$log2FC>=0, "enriched", "depleted")
  
  
  gene_enrichment<-droplevels(gene_enrichment)
  
  highlightedGene <- filter(gene_enrichment, gene == "kuz")
  highlightedGene <- droplevels(highlightedGene)
  
  p<-ggplot(gene_enrichment)
  p<-p + geom_bar(aes(gene, log2FC, fill = as.character(test)), stat="identity")
  #p<-p + geom_bar(data=highlightedGene, aes(gene, log2FC, fill="red"), colour="black", stat="identity")
  p<-p + guides(fill=FALSE)
  p<-p + scale_x_discrete("Gene")
  p<-p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 90, hjust=1),
          axis.text = element_text(size=7)
    )
  #p<-p + coord_flip()
  #p<-p + scale_y_reverse()
  
  gene_enrichment_plot <- paste("gene_enrichment.pdf")
  cat("Writing file", gene_enrichment_plot, "\n")
  ggsave(paste("plots/", gene_enrichment_plot, sep=""), width = 25, height = 5)
  
  
}



geneLenPlot <- function(n=0,gene_lengths_in="data/gene_lengths.txt"){
  gene_enrichment<-geneEnrichment(n=n)
  gene_lengths<-read.delim(gene_lengths_in, header = T)

  gene_enrichment$length<-as.numeric(gene_enrichment$length)
  gene_enrichment$log2FC<-as.numeric(gene_enrichment$log2FC)
  gene_enrichment$fc<-as.numeric(gene_enrichment$fc)
  gene_enrichment$observed<-as.numeric(gene_enrichment$observed)
  
  
  #  Var in x explained by Y
  # par(mfrow=c(2,2))
  # plot(enrichment_lm)
  
  # Set new col 'col' to indicate enrichment/depletion
  gene_enrichment$col<-as.factor(ifelse(gene_enrichment$log2FC > 0, 'enrichment', 'depletion'))
  
  # Only keep relevant cols
  gene_enrichment<-gene_enrichment[,c("gene","fpkm","observed",'expected','fc', 'log2FC', 'col')]
  gene_enrichment<-droplevels(gene_enrichment)
  
  # Join both df on 'gene'
  gene_lengths_df<-join(as.data.frame(gene_lengths), gene_enrichment, 'gene', type = "left")
  
  # Clean up null/na vals
  gene_lengths_df$fpkm <- ifelse(gene_lengths_df$fpkm=='NULL' | gene_lengths_df$fpkm=='NA' | is.na(gene_lengths_df$fpkm), 0, gene_lengths_df$fpkm)
  gene_lengths_df$level <- ifelse(gene_lengths_df$fpkm == 0 , 'not_expressed', 'expressed')
  gene_lengths_df$observed <- ifelse(gene_lengths_df$observed == 'NULL' | gene_lengths_df$observed == 'NA', 0, gene_lengths_df$observed)
  gene_lengths_df$level<-as.factor(gene_lengths_df$level)
  
  # Allow colouring of expressed/enriched/depleted 
  gene_lengths_df$col <- ifelse(is.na(gene_lengths_df$col), 'NA', as.character(gene_lengths_df$col))
  gene_lengths_df$col <- ifelse(gene_lengths_df$fpkm > 0, 'expressed', gene_lengths_df$col)
  
  # New col log10 length
  gene_lengths_df$log10length<-log10(gene_lengths$length)
  
  gene_lengths_df<-filter(gene_lengths_df, length >= 1000, length < 200000)
  gene_lengths_df<-droplevels(gene_lengths_df)
  # 
  # Linear model (predicter ~ predictor)
  enrichment_lm <- lm(observed ~ length, data = gene_lengths_df)
  # Exponential model
  enrichment_exp <- lm(observed^2 ~ length, data = gene_lengths_df)
  
  # plot( observed ~ length, data = gene_enrichment)
  lmRsq<-round(summary(enrichment_lm)$adj.r.squared, 2)
  expRsq<-round(summary(enrichment_exp)$adj.r.squared, 2)
  
  summary(pois <- glm(observed ~ length, family="poisson", data=gene_lengths_df))
  
  gene_lengths_p<-filter(gene_lengths_df, col != 'NA') 
  
  p <- ggplot(gene_lengths_p, aes(log10length, observed))
  p <- p + geom_jitter(aes(size=abs(log2FC), colour = col, alpha = 0.8))
  p <- p + scale_color_manual(values=c("#F8766D", "#00AFBB", "#E7B800"))
  p <- p + scale_x_continuous("Log10 Kb", limits=c(3, max(gene_lengths_p$log10length)))
  p <- p + scale_y_continuous("Count", limits=c(0,max(gene_lengths_p$observed)))
  p <- p + scale_size_continuous(range=c(0, abs(max(gene_lengths_p$log2FC))))

  p <- p + annotate(x = 3.5, y = 20, geom="text", label = paste('Lin:R^2:', lmRsq), size = 7,parse = TRUE)
  #p <- p + annotate(x = 30000, y = 18, geom="text", label = paste('Exp:R^2:', expRsq), size = 7,parse = TRUE)
  
  # Default model is formula = y ~ x
  # How much variation in X is explained by Y
  # How muc var in length is explained by observation
  
  p <- p + geom_smooth(method=lm, show.legend = FALSE) # linear
  #p <- p + geom_smooth(method=lm, formula = y ~ poly(x, 2), colour = "orange", show.legend = FALSE) #Quadratic
  #p <- p + geom_smooth(method=glm, method.args = list(family = "poisson"), colour = "red", se=T)
  p <- p + geom_smooth(colour="orange") # GAM
  p <- p + cleanTheme()
  p <- p + geom_rug(aes(colour=col,alpha=.8),sides="b")
  
  p <- p + guides(alpha = FALSE)
  
  colours<-c( "#E7B800", "#00AFBB")
  p2<-ggplot(gene_lengths_df)
  p2<-p2 + geom_density(aes(log10length, fill=level),alpha = 0.4)
  p2 <- p2 + cleanTheme()
  p2 <- p2 + scale_x_continuous("Log10 Kb", limits=c(3, max(gene_lengths_df$log10length)))
  p2 <- p2 + guides(alpha = FALSE)
  #p2 <- p2 + geom_rug(inherit.aes = F, aes(log10length,colour=level),alpha=0.2, sides = "tb")
  
  
  p2 <- p2 + geom_rug(data=subset(gene_lengths_df,level=="expressed"), aes(log10length,colour=level),alpha=0.7, sides = "b")
  p2 <- p2 + geom_rug(data=subset(gene_lengths_df,level=="not_expressed"), aes(log10length,colour=level),alpha=0.2, sides = "t")
  
  
  p2 <- p2 + scale_fill_manual(values=colours)
  p2 <- p2 + scale_colour_manual(values=colours)
  
  # p<-ggscatter(gene_lengths_p, x = "log10length", y = "observed",
  #           color = "col", size = "log2FC"
  #           )

  # p2<-ggdensity(gene_lengths_df, x = "log10length",
  #           add = "mean", rug = TRUE,
  #           color = "level", fill = "level",
  #           palette = c("#00AFBB", "#E7B800")
  #           )
  # 
  combined_plots <- ggarrange(p, p2, 
            labels = c("A", "B"),
            ncol = 1, nrow = 2)
  
  gene_len<-paste("gene_lengths_count_model_log10.pdf")
  cat("Writing file", gene_len, "\n")
  ggsave(paste("plots/", gene_len, sep=""), width = 10, height = 10)
  
  combined_plots
}

#' snvinGene
#'
#' Plot all snvs found in a given gene
#' @description Plot all snvs found in a given gene.
#' A 'gene_lengths' file must be provided with the following fields (cols 1..6 required)
#' gene length chrom    start      end      tss scaling_factor
#' This can be generated using the script 'script/genomic_features.pl' and a genome .gtf file
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
  wChrom<-as.character(region$chrom)
  wTss<-suppressWarnings(as.numeric(levels(region$tss))[region$tss])
  snv_data<-getData()
  snv_data<-filter(snv_data, chrom == wChrom & pos >= wStart & pos <= wEnd)
  
  if(nrow(snv_data) == 0){
    stop(paste("There are no snvs in", gene2plot, "- Exiting", "\n"))
  }
  
  p<-ggplot(snv_data)
  p<-p + geom_point(aes(pos/1000000, sample, colour = trans, size = 1.5), position=position_jitter(width=0, height=0.05))
  p<-p + guides(size = FALSE, sample = FALSE)
  p<-p + cleanTheme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.y = element_text(size = 30)
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(round(wStart/1000000, digits = 2),round(wEnd/1000000, digits = 2),by=0.05), limits=c(wStart/1000000, wEnd/1000000))
  p<-p + annotate("rect", xmin=region$start/1000000, xmax=region$end/1000000, ymin=0, ymax=0.3, alpha=.2, fill="skyblue")
  p<-p + geom_vline(xintercept = wTss/1000000, colour="red", alpha=.7, linetype="solid")
  
  p<-p + geom_segment(aes(x = wTss/1000000, y = 0, xend= wTss/1000000, yend = 0.1), colour="red")
  middle<-((wEnd/1000000+wStart/1000000)/2)
  p <- p + annotate("text", x = middle, y = 0.15, label=gene2plot, size=6)
  p<-p + ggtitle(paste("Chromosome:", wChrom))
  
  p
}


#' featuresHit
#'
#' Show top hit features
#' @import ggplot2
#' @keywords features
#' @export

featuresHit <- function(){
  snv_data<-getData()
  
  # To condense exon counts into "exon"
  snv_data$feature<-as.factor(gsub("exon_.*", "exon", snv_data$feature))
  
  # Reoders descending
  snv_data$feature<-factor(snv_data$feature, levels = names(sort(table(snv_data$feature), decreasing = TRUE)))
  
  #cols<-setCols(snv_data, "feature")
  
  p<-ggplot(snv_data)
  p<-p + geom_bar(aes(feature, fill = feature))
  #p<-p + cols
  p<-p + cleanTheme() +
    theme(axis.title.x=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous(expand = c(0.01, 0.01))
  
  # colour to a pub palette:
  # p<-p + ggpar(p, palette = 'jco')
  
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
  snv_data<-getData()
  snv_data<-filter(snv_data, gene != "intergenic")
  
  hit_count<-as.data.frame(sort(table(unlist(snv_data$gene)), decreasing = T))
  
  colnames(hit_count)<- c("gene", "count")
  head(hit_count, n)
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
  snv_data<-as.data.frame(bsapply(params))
  snv_data$genome<-as.integer(rowSums(snv_data))
  snv_data$genome_adj<-(snv_data$genome*2)
  
  if(!is.na(count)){
    tri_count<-snv_data['genome_adj']
    tri_count<-cbind(tri = rownames(tri_count), tri_count)
    colnames(tri_count) <- c("tri", "count")
    rownames(tri_count) <- NULL
    return(tri_count)
  }
  else{
    snv_data$x <- (1/snv_data$genome)
    scaling_factor<-snv_data['x']
    return(scaling_factor)
  }
  
}




#' tssDist
#'
#' Plot distance to TSS distribution
#' @param tss_pos File containing all TSS positions in genome: "gene chrom tss" [Default 'data/tss_positions.txt']
#' @param sim Simulate random SNVs accross genomic intervals? [Default: NO]
#' @param print Write the simulated random SNVs to a bed file ('data/simulatedSNVs.bed')? [Default: NO]
#' @import ggplot2
#' @keywords tss
#' @export

# tssDist <- function(tss_pos="data/tss_positions.txt",sim=NA, print=0){
#   tss_locations<-read.delim(tss_pos, header = T)
#   tss_locations$tss<-as.integer(tss_locations$tss)
#   snv_data<-getData()
#   if(!is.na(sim)){
#     simrep<-nrow(snv_data)
#     cat("Generating simulated snv_data for", simrep, "SNVs", "\n")
#     snv_data<-snvSim(N=simrep, write=print)
#     colnames(snv_data)<-c("chrom", "pos", "v3", "v4", "v5")
#     snv_data<-filter(snv_data, chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" | chrom == "Y" | chrom == "4")
#     snv_data<-droplevels(snv_data)
#   }
#   
#   fun2 <- function(p) {
#     index<-which.min(abs(tss_df$tss - p))
#     closestTss<-tss_df$tss[index]
#     chrom<-as.character(tss_df$chrom[index])
#     gene<-as.character(tss_df$gene[index])
#     dist<-(p-closestTss)
#     list(p, closestTss, dist, chrom, gene)
#   }
#   
#   l <- list()
#   
#   for (c in levels(snv_data$chrom)){
#     df<-filter(snv_data, chrom == c)
#     tss_df<-filter(tss_locations, chrom == c)
#     dist2tss<-lapply(df$pos, fun2)
#     dist2tss<-do.call(rbind, dist2tss)
#     dist2tss<-as.data.frame(dist2tss)
#     
#     colnames(dist2tss)=c("snp", "closest_tss", "min_dist", "chrom", "closest_gene")
#     dist2tss$min_dist<-as.numeric(dist2tss$min_dist)
#     l[[c]] <- dist2tss
#   }
#   
#   dist2tss<-do.call(rbind, l)
#   dist2tss<-as.data.frame(dist2tss)
#   dist2tss$chrom<-as.character(dist2tss$chrom)
#   
#   dist2tss<-arrange(dist2tss,(abs(min_dist)))
#   
#   # Removes chroms with fewer than 20 observations
#   snvCount <- table(dist2tss$chrom)
#   dist2tss <- subset(dist2tss, chrom %in% names(snvCount[snvCount > 25]))
#   dist2tss <- filter(dist2tss, chrom != 'Y')
#   
#   p<-ggplot(dist2tss)
#   p<-p + geom_density(aes(min_dist, fill = chrom), alpha = 0.3)
#   p<-p + scale_x_continuous("Distance to TSS (Kb)",
#                             limits=c(-100000, 100000),
#                             breaks=c( -10000, -1000, 0, 1000, 10000 ),
#                             expand = c(.0005, .0005),
#                             labels=c("-10", "-1", 0, "1", "10") )
#   p<-p + scale_y_continuous("Density", expand = c(0, 0))
#   p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
#   p<-p + cleanTheme()
#   #p<-p + facet_wrap(~chrom, ncol = 2, scales = "free_y")
#   p
#   
#   # p<-ggplot(dist2tss)
#   # p<-p + geom_histogram(aes(min_dist, fill = chrom), alpha = 0.6, bins=500)
#   # p<-p + scale_x_continuous("Distance to TSS", limits=c(-1000, 1000))
#   # p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
#   # p
#   
# }



tssDist <- function(tss_pos="data/tss_positions.txt",sim=NA, print=0,return=0){
  tss_locations<-read.delim(tss_pos, header = T)
  tss_locations$tss<-as.integer(tss_locations$tss)
  
  if(is.na(sim)){
    snv_data<-getData()
  }
  
  else{
    cat("Generating simulated snv_data\n")
    hit_count<-nrow(getData())
    snv_data<-snvSim(N=hit_count, write=print)
    colnames(snv_data)<-c("chrom", "pos", "v3", "v4", "v5")
    snv_data<-filter(snv_data, chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" | chrom == "Y" | chrom == "4")
    snv_data<-droplevels(snv_data)
  }
  
  
  
  # Will throw error if SVs don't exist on a chrom...
  
  # Removes chroms with fewer than 20 observations
  svCount <- table(snv_data$chrom)
  snv_data <- subset(snv_data, chrom %in% names(svCount[svCount > 30]))
  snv_data<-droplevels(snv_data)
  
  tss_locations <- subset(tss_locations, chrom %in% levels(snv_data$chrom))
  tss_locations<-droplevels(tss_locations)  
  
  
  fun2 <- function(p) {
    index<-which.min(abs(tss_df$tss - p))
    closestTss<-tss_df$tss[index]
    chrom<-as.character(tss_df$chrom[index])
    gene<-as.character(tss_df$gene[index])
    dist<-(p-closestTss)
    list(p, closestTss, dist, chrom, gene)
  }
  
  l <- list()
  
  for (c in levels(snv_data$chrom)){
    df<-filter(snv_data, chrom == c)
    tss_df<-filter(tss_locations, chrom == c)
    dist2tss<-lapply(df$pos, fun2)
    dist2tss<-do.call(rbind, dist2tss)
    dist2tss<-as.data.frame(dist2tss)
    
    colnames(dist2tss)=c("bp", "closest_tss", "min_dist", "chrom", "closest_gene")
    dist2tss$min_dist<-as.numeric(dist2tss$min_dist)
    l[[c]] <- dist2tss
  }
  
  dist2tss<-do.call(rbind, l)
  dist2tss<-as.data.frame(dist2tss)
  dist2tss$chrom<-as.character(dist2tss$chrom)
  
  dist2tss<-arrange(dist2tss,(abs(min_dist)))
  
  # Removes chroms with fewer than 20 observations
  # svCount <- table(dist2tss$chrom)
  # dist2tss <- subset(dist2tss, chrom %in% names(svCount[svCount > 10]))
  
  if(return==1){
    return(dist2tss)
  }
  else{
    p<-ggplot(dist2tss)
    p<-p + geom_density(aes(min_dist, fill = chrom), alpha = 0.3)
    p<-p + scale_x_continuous("Distance to TSS (Kb)",
                              limits=c(-10000, 10000),
                              breaks=c(-10000,-1000, 1000, 10000),
                              expand = c(.0005, .0005),
                              labels=c("-10", "-1", "1", "10") )
    p<-p + scale_y_continuous("Density")
    p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
    #p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 5)
    p <- p + geom_rug(aes(min_dist, colour=chrom))
    p<-p + cleanTheme() +
      theme(strip.text = element_text(size=20),
            legend.position="top")
    p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
    
    if(is.na(sim)){
      tssDistout<-paste("bpTSSdist.pdf")
    }
    else{
      tssDistout<-paste("bpTSSdist_sim.pdf")
    }
    cat("Writing file", tssDistout, "\n")
    ggsave(paste("plots/", tssDistout, sep=""), width = 20, height = 10)
    
    p
  }
}


tssDistOverlay <- function(){
  real_data<-tssDist(return=1)
  real_data$Source<-"Real"
  sim_data<-bpTssDist(sim=1, return=1)
  sim_data$Source<-"Sim"
  
  sim_data<-filter(sim_data, chrom != "Y", chrom != 4)
  sim_data<-droplevels(sim_data)
  real_data<-filter(real_data, chrom != "Y", chrom != 4)
  real_data<-droplevels(real_data)
  
  
  colours<-c( "#E7B800", "#00AFBB")
  
  
  p<-ggplot()
  p<-p + geom_density(data=real_data,aes(min_dist, fill = Source), alpha = 0.4)
  p<-p + geom_density(data=sim_data,aes(min_dist, fill = Source), alpha = 0.4)
  p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
  
  p<-p + scale_x_continuous("Distance to TSS (Kb)",
                            limits=c(-100000, 100000),
                            breaks=c(-100000,-10000,-1000, 1000, 10000, 100000),
                            expand = c(.0005, .0005),
                            labels=c("-100", "-10", "-1", "1", "10", "100") )
  p<-p + scale_y_continuous("Density")
  p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
  #p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 5)
  
  p <- p + geom_rug(data=real_data,aes(min_dist, colour=Source),sides="b")
  p <- p + geom_rug(data=sim_data,aes(min_dist, colour=Source),sides="t")
  
  
  p <- p + scale_fill_manual(values=colours)
  p <- p + scale_colour_manual(values=colours)
  
  p<-p + cleanTheme() +
    theme(strip.text = element_text(size=20),
          legend.position="top")
  
  p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
  
  
  tssDistout<-paste("bpTSSdist_overlay.pdf")
  
  cat("Writing file", tssDistout, "\n")
  ggsave(paste("plots/", tssDistout, sep=""), width = 20, height = 10)
  
  
  
  p
  
  
}


#' snvSim
#'
#' Generate simulated SNV hits acroos genomic regions (e.g. mappable regions)
#' @param intervals File containing genomic regions within which to simulate SNVs [Default 'data/intervals.bed]
#' @param N Number of random SNVs to generate [Default nrow(snv_data)]
#' @import GenomicRanges
#' @keywords sim
#' @export

snvSim <- function(intervals="data/intervals.bed", N=1000, write=F){
  suppressPackageStartupMessages(require(GenomicRanges))
  
  intFile <- import.bed(intervals)
  space <- sum(width(intFile))
  positions <- sample(c(1:space), N)
  new_b <- GRanges(seqnames=as.character(rep(seqnames(intFile), width(intFile))),
                   ranges=IRanges(start=unlist(mapply(seq, from=start(intFile), to=end(intFile))), width=1))
  bedOut<-new_b[positions]
  if(write){
    export.bed(new_b[positions], "data/simulatedSNVs.bed")
  }
  remove(new_b)
  return(data.frame(bedOut))
}



svDist <- function(svs="data/all_bps_new.txt",sim=NA, print=0){
  svBreaks<-read.delim(svs, header = F)
  colnames(svBreaks) <- c("event", "bp_no", "sample", "chrom", "bp", "gene", "feature", "type", "length")
  
  svBreaks$bp<-as.integer(svBreaks$bp)
  svBreaks<-filter(svBreaks, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  svBreaks <- droplevels(svBreaks)
  
  snv_data<-getData()
  
  if(!is.na(sim)){
    simrep<-nrow(snv_data)
    cat("Generating simulated snv_data for", simrep, "SNVs", "\n")
    snv_data<-snvSim(N=simrep, write=print)
    snv_data$end<-NULL
    snv_data$width<-NULL
    snv_data$strand<-NULL
    snv_data$sample <- as.factor(sample(levels(svBreaks$sample), size = nrow(snv_data), replace = TRUE))
    #snv_data$type <- sample(levels(svBreaks$type), size = nrow(snv_data), replace = TRUE)
    colnames(snv_data)<-c("chrom", "pos", "sample")
    snv_data<-filter(snv_data, chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" | chrom == "Y" | chrom == "4")
    snv_data<-droplevels(snv_data)
  }
  
  snv_data <- subset(snv_data, sample %in% levels(svBreaks$sample))
  snv_data <- droplevels(snv_data)
  
  snv_data <- subset(snv_data, chrom %in% levels(svBreaks$chrom))
  snv_data <- droplevels(snv_data)
  
  
  fun3 <- function(p) {
    index<-which.min(abs(sv_df$bp - p))
    closestBp<-as.numeric(sv_df$bp[index])
    chrom<-as.character(sv_df$chrom[index])
    gene<-as.character(sv_df$gene[index])
    sample<-as.character(sv_df$sample[index])
    type<-as.character(sv_df$type[index])
    
    dist<-(p-closestBp)
    list(p, closestBp, dist, chrom, gene, type, sample)
  }
  
  l <- list()
  
  for (c in levels(snv_data$chrom)){
    for (s in levels(snv_data$sample)){
      df<-filter(snv_data, chrom == c & sample == s)
      sv_df<-filter(svBreaks, chrom == c & sample == s)
      
      dist2bp<-lapply(df$pos, fun3)
      dist2bp<-do.call(rbind, dist2bp)
      dist2bp<-as.data.frame(dist2bp)
      
      colnames(dist2bp)=c("snp", "closest_bp", "min_dist", "chrom", "closest_gene", "type", "sample")
      dist2bp$min_dist<-as.numeric(dist2bp$min_dist)
      l[[s]] <- dist2bp
    }
    l[[c]] <- dist2bp
  }
  
  dist2bp<-do.call(rbind, l)
  dist2bp<-as.data.frame(dist2bp)
  dist2bp$chrom<-as.character(dist2bp$chrom)
  dist2bp$type<-as.character(dist2bp$type)
  
  snvCount <- table(dist2bp$chrom)
  dist2bp <- subset(dist2bp, chrom %in% names(snvCount[snvCount > 25]))
  
  dist2bp<-arrange(dist2bp,(abs(min_dist)))
  
  dist2bp <- dist2bp %>% na.omit()
  
  p<-ggplot(dist2bp)
  p<-p + geom_density(aes(min_dist, fill=type), alpha = 0.3)
  # p<-p + scale_x_continuous("Distance to SV BP (Kb)",
  #                           limits=c(-10000000, 10000000),
  #                           breaks=c(-10000000, -100000, -10000, -1000, 0, 1000, 10000, 100000, 10000000),
  #                           expand = c(.0005, .0005),
  #                           labels=c("-10000", "-100", "-10", "-1", 0, "1", "10", "100", "10000") )
  p<-p + scale_y_continuous("Density", expand = c(0, 0))
  p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted") 
  p<-p + cleanTheme()
  p<-p + facet_wrap(~chrom, ncol = 2, scales = "free_y")
  p
  
  # p<-ggplot(dist2bp)
  # p<-p + geom_histogram(aes(as.numeric(min_dist, fill = type)), alpha = 0.6, binwidth = 1000)
  # p<-p + scale_x_continuous("Distance to TSS", limits=c(-1000000, 1000000))
  # p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
  # p
  
}



#' chromDist
#'
#' Plot genome-wide snv distribution 
#' @import ggplot2
#' @keywords distribution
#' @export

chromDist <- function(object=NA, notch=0){
  snv_data<-getData()
  ext<-'.pdf'
  if(is.na(object)){
    object<-'grouped_trans'
    cols<-setCols(snv_data, "grouped_trans")
  }
  
  if(notch){
    snv_data<-exclude_notch()
    ext<-'_excl.N.pdf'
  }
  
  cat("Plotting snvs by", object, "\n")
  
  p<-ggplot(snv_data)
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


###########
##  Misc ##
###########

# Gene lengths on accross chroms

#gene_lengths="data/gene_lengths.txt"
#gene_lengths<-read.delim(gene_lengths, header = T)

#p<-ggplot(gene_lengths)
#p<-p + geom_point(aes( x=(start+(length/2)), y=log10(length), colour = chrom))
#p<-p + facet_wrap(~chrom, scale = "free_x")
#p

