#' Merge old and new records
#'
#' Merge the old and new (re-drilled) oxygen isotope records of a specific site,
#' either at the start or the end of the overlap period and with or without
#' adjusting the mean values of the records before merging.
#'
#' @param site character string of the site ID for which to merge the old and
#'   new records; must match the list of sites for which re-drilled records are
#'   available.
#' @param data data frame with the NGT and associated oxygen isotope records.
#' @param method integer signalling the merging method: for \code{method = 1}
#'   the new record is used over the entire length of the common overlap period
#'   until its start (year CE), for \code{method = 2} it is used only until the
#'   end of the common overlap period.
#' @param adjustMean logical; shall the mean values of the two records within
#'   the common overlap period be adjusted prior to the merging? Defaults to
#'   \code{FALSE}. If \code{TRUE}, mean adjustment is always performed for the
#'   new record relative to the old one.
#' @return numeric vector with the merged record.
#' @author Thomas Münch
doMerge <- function(site = "B18", data, method = 1, adjustMean = FALSE) {

  if (!site %in% c("B18", "B21", "B23", "B26", "NGRIP", "stack"))
    stop("Unsuitable site requested.")

  if (!method %in% c(1, 2))
    stop("Only 'method = 1' or 'method = 2' available.")

  startYear <- data$Year[1]
  cores <- calculateOverlapStatistics(data = data, site = site)

  if (adjustMean) {  
    cores$dat$pairData[, 3] <-
      cores$dat$pairData[, 3] - cores$stat$diffMeanOverlap
  }

  if (method == 1) {

    i <- match(startYear : cores$stat$startOverlap, data$Year)
    merged <- c(cores$dat$pairData[i, 3], cores$dat$pairData[-i, 2])

  } else {

    i <- match(startYear : (cores$stat$endOverlap + 1), data$Year)
    merged <- c(cores$dat$pairData[i, 3], cores$dat$pairData[-i, 2])

  }

  return(merged)

}

#' Merge a set of old and new records
#'
#' Wrapper function to merge all suitable pairs of old and new oxygen isotope
#' records from the NGT data set.
#'
#' @param data data frame with the NGT and associated oxygen isotope records.
#' @param sites vector of the character IDs of the sites for which the merging
#'   shall be conducted; defaults to the list of sites with available paired
#'   records.
#' @param method integer signalling the merging method: for \code{method = 1}
#'   the new record is used over the entire length of the common overlap period
#'   until its start (year CE), for \code{method = 2} it is used only until the
#'   end of the common overlap period.
#' @param adjustMean logical; shall the mean values of the two records within
#'   the common overlap period be adjusted prior to the merging? Defaults to
#'   \code{FALSE}. If \code{TRUE}, mean adjustment is always performed for the
#'   new record relative to the old one.
#' @return a data frame with an age column and the merged versions of the
#'   records specified in \code{sites}.
#' @author Thomas Münch
mergeCores <- function(data,
                       sites = c("B18", "B21", "B23", "B26", "NGRIP", "stack"),
                       method = 1, adjustMean = FALSE) {

  if (!method %in% c(1, 2))
    stop("Only 'method = 1' or 'method = 2' available.")

  if (any(is.na(match(sites, names(data)))))
    stop("Requested site(s) not found in given data.")

  merged <- sapply(sites, doMerge, data, method, adjustMean)
  merged <- as.data.frame(cbind(Year = data$Year, merged))

  return(merged)

}
