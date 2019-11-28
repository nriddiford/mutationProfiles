#' indel_count
#'
#' Plot indel size per sample
#' @import dplyr
#' @keywords count
#' @export
#'
plot_indels <- function(..., indel_data=NULL){

  if(missing(indel_data)){
    indel_data <- mutationProfiles::getData(..., infile ='data/annotated_indels.txt', exclude=F, type = 'indel')
  }

  calls_by_sample <- indel_data %>%
    dplyr::group_by(sample) %>%
    dplyr::tally()

  blueBar <- '#3B8FC7'
  calls_by_sample$colour <- blueBar

  calls_by_sample <- transform(calls_by_sample, sample = reorder(sample, -n))

  p <- ggplot(calls_by_sample, aes(sample, n, fill = colour))
  p <- p + geom_bar(stat='identity')
  p <- p + scale_y_continuous("Number of calls", limits=c(0,max(calls_by_sample$n)), breaks=seq(0,max(calls_by_sample$n), by=100))
  p <- p + scale_x_discrete("Sample")
  p <- p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
          axis.text = element_text(size=20), axis.title = element_text(size=20)
    )
  p <- p + scale_fill_identity()
  p
}

#' indel_length
#'
#' Plot indel size per sample
#' @import dplyr
#' @keywords count
#' @export
indel_lengths <- function(..., indel_data=NULL, plot=TRUE, count=TRUE){

  if(missing(indel_data)){
    indel_data <- mutationProfiles::getData(..., infile ='data/annotated_indels.txt', exclude=F, type = 'indel')
  }
  indel_lengths <- indel_data %>%
    dplyr::mutate(length = nchar(as.character(alt))) %>%
    dplyr::mutate(type_class = as.factor(ifelse(length>3,
                                      paste0(tolower(type), '>', '3bp'),
                                      paste0(tolower(type), '_', length, 'bp')
                                      ))) %>%
    dplyr::group_by(sample, type_class) %>%
    dplyr::tally()


  indel_lengths <- indel_lengths %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(total = sum(n)) %>%
    dplyr::mutate(fraction = n/total)

  if(plot){
    if(count) {
      p <- ggplot(indel_lengths, aes(fct_reorder(sample, -n), n, fill=type_class))
      ytitle <- "Number of calls"
    } else {
      p <- ggplot(indel_lengths, aes(fct_reorder(sample, -n), fraction, fill=type_class))
      ytitle <- "Fraction of calls"
    }

    p <- p + geom_bar(stat='identity')
    p <- p + scale_y_continuous("Number of calls", expand = c(0,0))
    p <- p + scale_x_discrete("Sample")
    p <- p + cleanTheme() +
      theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
            axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
            axis.text = element_text(size=20), axis.title = element_text(size=20)
      )
    p
  } else {
    return(indel_lengths)
  }
}


#' indel_seqs
#'
#' Tally inserted/deleted sequences
#' @import dplyr
#' @export
indel_seqs <- function(..., indel_data=NULL, plot=TRUE){

  if(missing(indel_data)){
    indel_data <- mutationProfiles::getData(..., infile ='data/annotated_indels.txt', exclude=F, type = 'indel')
  }
  seqs <- indel_data %>%
    dplyr::mutate(sex = ifelse(stringr::str_detect(sample, 'D'), "F", "M")) %>%
    dplyr::group_by(sex, type, alt) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(freq = (n / sum(n)) * 100) %>%
    dplyr::filter(freq > 1)

  if(plot){
    p <- ggplot(seqs, aes(fct_reorder(alt, freq), freq, fill=type))
    p <- p + geom_bar(stat='identity')
    p <- p + scale_y_continuous("Percentage", expand = c(0,0))
    p <- p + scale_x_discrete("Sequence")
    p <- p + cleanTheme() +
      theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
            axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
            axis.text = element_text(size=20), axis.title = element_text(size=20)
      )
    p <- p + facet_wrap(sex~type)
    p <- p + coord_flip()
    p
  } else {
    return(indel_lengths)
  }
}
