library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(deconstructSigs)
library(reshape)
library(data.table)


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
  
  colnames(data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "a_freq", "caller", "feature", "gene")
  data$a_freq<-as.numeric(levels(data$a_freq))[data$a_freq]
  
  #data<-filter(data, is.na(a_freq))
  #data<-filter(data, a_freq >= 0.20)

  #filter out samples
  data<-filter(data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  #data<-filter(data, sample != "A373R11" & sample != 'A373R13')
  
  
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
#' Get colours for n levels
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


#' rainfall
#'
#' Plot log10 distances between snvs as rainfall plot
#' @import ggplot2
#' @keywords rainfall
#' @export

rainfall <- function(){
  data<-getData()

  data<-filter(data, sample != "A373R11" & sample != 'A373R13')
  distances<-do.call(rbind, lapply(split(data[order(data$chrom, data$pos),], data$chrom[order(data$chrom, data$pos)]),
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
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
          )
  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  
  rainfall_out<-paste("rainfall.pdf")
  cat("Writing file", rainfall_out, "\n")
  ggsave(paste("plots/", rainfall_out, sep=""), width = 20, height = 10)
  
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
  
  cols <- setCols(mutData, 'signature')
  
  p <- ggplot(mutData[order(mutData$signature),])
  p <- p + geom_bar(aes(reorder(sample, -score), score, fill=signature),colour="black", stat = "identity")
  p <- p + scale_x_discrete("Sample")
  p <- p + scale_y_continuous("Signature contribution", expand = c(0.01, 0.01), breaks=seq(0, 1, by=0.1))
  p <- p + cleanTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  p <- p + cols
  p
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
  data<-getData()
  data<-filter(data, chrom == wChrom & pos >= wStart & pos <= wEnd)
  
  if(nrow(data) == 0){
    stop(paste("There are no snvs in", gene2plot, "- Exiting", "\n"))
  }
  
  p<-ggplot(data)
  p<-p + geom_point(aes(pos/1000000, sample, colour = trans, size = 1.5), position=position_jitter(width=0, height=0.05))
  p<-p + guides(size = FALSE, sample = FALSE)
  p<-p + cleanTheme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(round(wStart/1000000, digits = 2),round(wEnd/1000000, digits = 2),by=0.05), limits=c(wStart/1000000, wEnd/1000000))
  p<-p + annotate("rect", xmin=region$start/1000000, xmax=region$end/1000000, ymin=0, ymax=0.1, alpha=.2, fill="skyblue")
  p<-p + geom_vline(xintercept = wTss/1000000, colour="red", alpha=.7, linetype="solid")
  
  p<-p + geom_segment(aes(x = wTss/1000000, y = 0, xend= wTss/1000000, yend = 0.1), colour="red")
  middle<-((wEnd/1000000+wStart/1000000)/2)
  p <- p + annotate("text", x = middle, y = 0.05, label=gene2plot, size=6)
  p<-p + ggtitle(paste("Chromosome:", wChrom))
  
  p
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

tssDist <- function(tss_pos="data/tss_positions.txt",sim=NA, print=0){
  tss_locations<-read.delim(tss_pos, header = T)
  tss_locations$tss<-as.integer(tss_locations$tss)
  data<-getData()
  if(!is.na(sim)){
    simrep<-nrow(data)
    cat("Generating simulated data for", simrep, "SNVs", "\n")
    data<-snvSim(N=simrep, write=print)
    colnames(data)<-c("chrom", "pos", "v3", "v4", "v5")
    data<-filter(data, chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" | chrom == "Y" | chrom == "4")
    data<-droplevels(data)
  }

  fun2 <- function(p) {
    index<-which.min(abs(tss_df$tss - p))
    closestTss<-tss_df$tss[index]
    chrom<-as.character(tss_df$chrom[index])
    gene<-as.character(tss_df$gene[index])
    dist<-(p-closestTss)
    list(p, closestTss, dist, chrom, gene)
  }
  
  l <- list()
  
  for (c in levels(data$chrom)){
    df<-filter(data, chrom == c)
    tss_df<-filter(tss_locations, chrom == c)
    dist2tss<-lapply(df$pos, fun2)
    dist2tss<-do.call(rbind, dist2tss)
    dist2tss<-as.data.frame(dist2tss)
    
    colnames(dist2tss)=c("snp", "closest_tss", "min_dist", "chrom", "closest_gene")
    dist2tss$min_dist<-as.numeric(dist2tss$min_dist)
    l[[c]] <- dist2tss
  }
  
  dist2tss<-do.call(rbind, l)
  dist2tss<-as.data.frame(dist2tss)
  dist2tss$chrom<-as.character(dist2tss$chrom)
  
  dist2tss<-arrange(dist2tss,(abs(min_dist)))
  
  # Removes chroms with fewer than 20 observations
  snvCount <- table(dist2tss$chrom)
  dist2tss <- subset(dist2tss, chrom %in% names(snvCount[snvCount > 25]))
  dist2tss <- filter(dist2tss, chrom != 'Y')
 
  p<-ggplot(dist2tss)
  p<-p + geom_density(aes(min_dist, fill = chrom), alpha = 0.3)
  p<-p + scale_x_continuous("Distance to TSS (Kb)",
                            limits=c(-100000, 100000),
                            breaks=c( -10000, -1000, 0, 1000, 10000 ),
                            expand = c(.0005, .0005),
                            labels=c("-10", "-1", 0, "1", "10") )
  p<-p + scale_y_continuous("Density", expand = c(0, 0))
  p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
  p<-p + cleanTheme()
  #p<-p + facet_wrap(~chrom, ncol = 2, scales = "free_y")
  p
  
  # p<-ggplot(dist2tss)
  # p<-p + geom_histogram(aes(min_dist, fill = chrom), alpha = 0.6, bins=500)
  # p<-p + scale_x_continuous("Distance to TSS", limits=c(-1000, 1000))
  # p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
  # p
  
}


#' snvSim
#'
#' Generate simulated SNV hits acroos genomic regions (e.g. mappable regions)
#' @param intervals File containing genomic regions within which to simulate SNVs [Default 'data/intervals.bed]
#' @param N Number of random SNVs to generate [Default nrow(data)]
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
  
  data<-getData()
  
  if(!is.na(sim)){
    simrep<-nrow(data)
    cat("Generating simulated data for", simrep, "SNVs", "\n")
    data<-snvSim(N=simrep, write=print)
    data$end<-NULL
    data$width<-NULL
    data$strand<-NULL
    data$sample <- as.factor(sample(levels(svBreaks$sample), size = nrow(data), replace = TRUE))
    #data$type <- sample(levels(svBreaks$type), size = nrow(data), replace = TRUE)
    colnames(data)<-c("chrom", "pos", "sample")
    data<-filter(data, chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" | chrom == "Y" | chrom == "4")
    data<-droplevels(data)
  }
  
  data <- subset(data, sample %in% levels(svBreaks$sample))
  data <- droplevels(data)
  
  data <- subset(data, chrom %in% levels(svBreaks$chrom))
  data <- droplevels(data)
  
  
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

  for (c in levels(data$chrom)){
    for (s in levels(data$sample)){
      df<-filter(data, chrom == c & sample == s)
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
#' @return A data frame with FC scores for all genes seen at least n times in snv data
#' @export 

featureEnrichment <- function(features='data/genomic_features.txt', genome_length=118274340){
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


EnrichmentPlot <- function() {
  feature_enrichment<-featureEnrichment()

  feature_enrichment$Log2FC <- log2(as.numeric(feature_enrichment$fc))
  
  feature_enrichment$feature <- as.character(feature_enrichment$feature)
  feature_enrichment$fc <- as.numeric(feature_enrichment$fc)
  
  feature_enrichment <- transform(feature_enrichment, feature = reorder(feature, -fc))
  
  feature_enrichment <- filter(feature_enrichment, observed >= 5)
  
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
#' @param n The number of times we need to have seen a gene in our data to view its enrichment score [Default 3]
#' @param genome_length The total legnth of the genome [Default 137547960 (chroms 2, 3, 4, X & Y for Drosophila melanogastor Dmel6.12)]
#' @keywords enrichment
#' @import dplyr
#' @return A data frame with FC scores for all genes seen at least n times in snv data
#' @export 

geneEnrichment <- function(gene_lengths="data/gene_lengths.txt", n=3, genome_length=118274340){
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
  genesFC<-filter(genesFC, observed >= n)
  # Sort by FC value
  genesFC<-arrange(genesFC,desc(as.integer(fc)))
  return(genesFC)
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


# svBreaks<-read.delim("../DUMMYSV.txt", header = F)
# colnames(svBreaks) <- c("event", "bp_no", "sample", "chrom", "bp", "gene", "feature", "type", "length")
# svBreaks<-head(svBreaks,10)
# svBreaks<-select(svBreaks, sample, chrom, bp, gene, type)
# svBreaks <- droplevels(svBreaks)

# 
# data<-read.delim("DUMMY_DATA.txt", header=F)
# colnames(data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "a_freq", "caller", "feature", "gene")
# data<-head(data,10)
# data<-select(data, sample, chrom, pos)
# 
# data <- subset(data, sample %in% levels(svBreaks$sample))
# data <- droplevels(data)
# 
# data <- subset(data, chrom %in% levels(svBreaks$chrom))
# data <- droplevels(data)

# svBreaks<-structure(list(sample = structure(c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 
#                                               1L, 2L, 1L), .Label = c("S1", "S2"), class = "factor"), chrom = structure(c(1L, 
#                                                                                                                           1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L), .Label = c("2L", "2R"), class = "factor"), 
#                          bp = c(2425901L, 2426025L, 6694426L, 6694566L, 8387755L, 
#                                 8387927L, 8963713L, 963799L, 980364L, 980521L), gene = structure(c(3L, 
#                                                                                                    3L, 5L, 5L, 4L, 4L, 2L, 2L, 1L, 1L), .Label = c("CG8213", 
#                                                                                                                                                    "CG8216", "intergenic", "pdm3", "Tsp"), class = "factor"), 
#                          type = structure(c(2L, 1L, 2L, 1L, 3L, 3L, 3L, 4L, 4L, 3L
#                          ), .Label = c("DEL", "DUP", "INV", "TANDUP"), class = "factor")), row.names = c(NA, 
#                                                                                                          10L), .Names = c("sample", "chrom", "bp", "gene", "type"), class = "data.frame")
# 
# 
# data<-structure(list(sample = structure(c(1L, 2L, 2L, 1L, 1L, 1L, 1L, 
#                                           2L, 2L, 2L), .Label = c("S1", "S2"), class = "factor"), chrom = structure(c(1L, 
#                                                                                                                       1L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 2L), .Label = c("2L", "2R"), class = "factor"), 
#                      pos = c(318351L, 605574L, 1014043L, 2031592L, 2886957L, 2910379L, 
#                              2218351L, 105574L, 1344043L, 216957L)), .Names = c("sample", 
#                                                                                 "chrom", "pos"), row.names = c(NA, 10L), class = "data.frame")
# 
# 
# 
# fun3 <- function(p) {
#   index<-which.min(abs(sv_df$bp - p))
#   closestBp<-as.numeric(sv_df$bp[index])
#   chrom<-as.character(sv_df$chrom[index])
#   gene<-as.character(sv_df$gene[index])
#   sample<-as.character(sv_df$sample[index])
#   type<-as.character(sv_df$type[index])
#   
#   dist<-(p-closestBp)
#   list(p, closestBp, dist, chrom, gene, type, sample)
# }
# 
# l <- list()
# 
# for (c in levels(data$chrom)){
#   for (s in levels(data$sample)){
# 
#     df<-filter(data, chrom == c & sample == s)
#     sv_df<-filter(svBreaks, chrom == c & sample == s)
#     
#     dist2bp<-lapply(df$pos, fun3)
#     dist2bp<-do.call(rbind, dist2bp)
#     dist2bp<-as.data.frame(dist2bp)
#     
#     colnames(dist2bp)=c("snp", "closest_bp", "min_dist", "chrom", "closest_gene", "type", "sample")
#     cat("Chrom = ", c, "Sample = ", s, "\n")
#     dist2bp$min_dist<-as.numeric(dist2bp$min_dist)
#     l[[s]] <- dist2bp
#   }
#   #l[[c]] <- dist2bp
# }
# 
# dist2bp<-do.call(rbind, l)
# dist2bp<-as.data.frame(dist2bp)
# dist2bp$chrom<-as.character(dist2bp$chrom)
# dist2bp$type<-as.character(dist2bp$type)
# 
# print(dist2bp)
# 


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

