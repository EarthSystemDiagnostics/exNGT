#' Calculate standard error
#'
#' Calculate the standard error of a data vector, i.e. its standard deviation
#' divided by the square root of the number of data points.
#'
#' @param x numeric vector for which to calculate the standard error.
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'   stripped before the computation proceeds.
#' @return the standard error of \code{x}.
#' @author Maria Hoerhold
#'
sdError <- function(x, na.rm = FALSE) {
  sd(x, na.rm = na.rm) / sqrt(length(x))
}

#' Statistics of the overlap period
#'
#' Deliver the statistics of the overlap period between the old and the
#' new (re-drilled) isotope record at a certain site.
#'
#' @param data data frame with the NGT and associated oxygen isotope records.
#' @param site character string of the site ID for which to calculate the
#'   overlap statistics; must match the list of sites for which re-drilled
#'   records are available.
#' @return a nested list with two elements:
#'   * "dat": a list of two data frames with the full data for the paired
#'     records at the requested site and the data for the overlap period.
#'   * "stat": a list of six elements with statistical information about the
#'     records' overlap; including start and end year of overlap period, mean
#'     isotope values in the overlap period, the difference of these mean
#'     values, the correlation for the overlap period, and the standard error of
#'     the mean values in the overlap period.
#' @author Thomas Münch, Maria Hoerhold
calculateOverlapStatistics <- function(data, site = "B18") {

  if (!site %in% c("B18", "B21", "B23", "B26", "NGRIP", "stack"))
    stop("Unsuitable site requested.")

  pairID <- c("Year", site, paste(site, "12", sep = "_"))

  pairData <- data[, pairID]
  overlap <- na.omit(pairData)
  meanOverlap <- colMeans(overlap[, -1])
  sdErrorOverlap <- apply(overlap[, -1], 2, sdError)
  
  output <- list(

    dat  = list(pairData = pairData, overlapData = overlap),
    stat = data.frame(
      corePair          = site,
      startOverlap      = min(overlap$Year),
      endOverlap        = max(overlap$Year),
      meanOverlapOld    = meanOverlap[1],
      meanOverlapNew    = meanOverlap[2],
      sdErrorOverlapOld = sdErrorOverlap[1],
      sdErrorOverlapNew = sdErrorOverlap[2],
      diffMeanOverlap   = diff(meanOverlap),
      corrOverlap       = cor(overlap[, 2], overlap[, 3]),
      row.names = NULL)
  )
  
  return(output)
  
}
