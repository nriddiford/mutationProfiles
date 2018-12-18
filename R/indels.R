#' indel_count
#'
#' Plot indel size per sample
#' @import dplyr
#' @keywords count
#' @export
#'
plot_indels <- function(..., indels='data/annotated_indels.txt'){

  indel_data <- mutationProfiles::getData(infile = indels, type = 'indel')

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
#'
plot_indel_lengths <- function(..., exclude=TRUE, indels='data/annotated_indels.txt'){

  indel_data <- mutationProfiles::getData(...,infile = indels, exclude=exclude, type = 'indel')

  # indel_data$length <- nchar(as.character(indel_data$alt))-3
  #
  # indel_lengths$type_class <- ifelse(indel_lengths$length>2,
  #                                    paste0(tolower(indel_lengths$type), '>', '2bp'),
  #                                    paste0(tolower(indel_lengths$type), '_', indel_lengths$length, 'bp')
  #                                    )

  indel_lengths <- indel_data %>%
    dplyr::mutate(length = nchar(as.character(alt))-3) %>%
    dplyr::mutate(type_class = as.factor(ifelse(length>2,
                                      paste0(tolower(type), '>', '2bp'),
                                      paste0(tolower(type), '_', length, 'bp')
                                      ))) %>%
    dplyr::group_by(sample, type_class) %>%
    dplyr::tally()

  indel_lengths <- indel_lengths %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(total = sum(n))

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
}
