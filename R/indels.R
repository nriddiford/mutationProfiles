#' indel_count
#'
#' Plot indel size per sample
#' @import dplyr
#' @keywords count
#' @export
#'
plot_indels <- function(..., exclude=TRUE, indels='data/annotated_indels.txt'){
  indel_data <- mutationProfiles::getData(..., infile = indels, exclude=exclude, type = 'indel')
  calls_by_sample <- indel_data %>%
    dplyr::group_by(sample) %>%
    dplyr::tally()

  calls_by_sample <- transform(calls_by_sample, sample = reorder(sample, -n))

  p <- ggplot(calls_by_sample, aes(sample, n))
  p <- p + geom_bar(stat='identity')
  p <- p + scale_y_continuous("Number of calls", limits=c(0,max(calls_by_sample$n)), breaks=seq(0,max(calls_by_sample$n), by=100))
  p <- p + scale_x_discrete("Sample")
  p <- p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
          axis.text = element_text(size=20), axis.title = element_text(size=20)
    )
  p
}

#' indel_length
#'
#' Plot indel size per sample
#' @import dplyr
#' @keywords count
#' @export
indel_lengths <- function(..., exclude=TRUE, indels='data/annotated_indels.txt', plot=TRUE){

  indel_data <- mutationProfiles::getData(...,infile = indels, exclude=exclude, type = 'indel')

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
    dplyr::mutate(total = sum(n))

  if(plot){
    p <- ggplot(indel_lengths, aes(fct_reorder(sample, -total), n, fill=type_class))
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
