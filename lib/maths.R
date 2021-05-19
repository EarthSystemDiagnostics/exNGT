#' Filter time series
#'
#' Apply a given filter to a time series using different endpoint constraints as
#' described below.
#'
#' Note that when passing objects of class \code{ts}, the time step provided is
#' not used; thus, for time series with a time step different from 1, the filter
#' has to be adapted accordingly.
#'
#' Leading and trailing NA values are automatically stripped from the input
#' vector so that they do not spread into the filtered data when applying the
#' endpoint constraints, but added in again after filtering so that the output
#' vector has the same length as the input. This does not apply to any internal
#' NA values, which instead are handled by \code{na.rm}.
#'
#' The function applies endpoint constrains following Mann et al., GRL, 2004;
#' available methods are:
#' \itemize{
#'   \item method = 0: no constraint (loss at both ends);
#'   \item method = 1: minimum norm constraint;
#'   \item method = 2: minimum slope constraint;
#'   \item method = 3: minimum roughness constraint;
#'   \item method = 4: circular filtering.
#' }
#'
#' @param data numeric vector with the input timeseries (standard or ts object).
#' @param filter numeric vector of filter weights.
#' @param method single integer for choosing an endpoint constraint method;
#'   available choices are integers 0-4.
#' @param na.rm logical; control the handling of internal NA values in
#'   \code{data}. If set to \code{TRUE}, any internal NA values are removed by
#'   linear interpolation from the neighbouring values; defaults to
#'   \code{FALSE}.
#' @return a ts object with the filtered timeseries.
#' @author Thomas Laepple
#' @source The endpoint constraint methods are based on the study:
#'   Michael E. Mann, On smoothing potentially non‐stationary climate time
#'   series, Geophys. Res. Lett., 31, L07214, doi:10.1029/2004GL019569, 2004.
#'
ApplyFilter <- function(data, filter, method = 0, na.rm = FALSE) {

  if (!method %in% (0 : 4))
    stop("Unknown method; only 0 : 4 available.")

  result <- rep(NA, length(data))

  # remove leading and trailing NA's
  x <- c(zoo::na.trim(data))
  n <- length(x)

  # linearly interpolate internal NA's if requested
  if (na.rm) {x <- stats::approx(1 : n, x, 1 : n)$y}

  circular = FALSE

  if (method == 0 | method == 4) {

    if (method == 4) {circular = TRUE}

    xf <- stats::filter(x, filter, circular = circular)

  } else {

    N <- floor(length(filter) / 2)

    if (method == 1) {

      before <- rep(mean(x), N)
      after  <- rep(mean(x), N)

    } else if (method == 2 | method == 3) {

      before <- x[N : 1]
      after  <- x[n : (n - N + 1)]

      if (method == 3) {

        before <- x[1] - (before - mean(before))
        after  <- x[n] - (after - mean(after))

      }
    }

    xf <- stats::filter(c(before, x, after), filter, circular = circular)
    xf <- xf[(N + 1) : (N + n)]

  }

  i <- seq(match(x[1], data), by = 1, length.out = n)
  result[i] <- xf

  return(ts(result, frequency = frequency(data)))

}

#' Calculate running mean of a set of time series
#'
#' Calculate the running mean of a set of time series, given as a data frame or
#' matrix, following different endpoint constraints.
#'
#' @param data a single numeric vector to filter with a running mean window, or
#'   a data frame or matrix for which the columns are filtered.
#' @param window single integer giving the size of the running mean window.
#' @param method single integer for choosing an endpoint constraint method;
#'   available choices are integers 0-4 (see \code{\link{ApplyFilter}} for
#'   details).
#' @param hasAgeColumn logical; do the supplied data contain an age column as
#'   the first column?
#' @return the running mean of the input \code{data}.
#' @author Maria Hoerhold
#'
filterData <-function(data, window = 5, method = 2, hasAgeColumn = TRUE) {

  if (window == 1) return(data)

  if (hasAgeColumn) {
    if (is.null(dim(data))) {
      stop("'data' is not a matrix or data frame.", call. = FALSE)
    }
    if (ncol(data) < 2) {
      stop("'data' only has 1 column; need at least 2.", call. = FALSE)
    }
    colNames <- colnames(data)
    Year <- data[, 1]
    data <- data[, -1]
  }
  
  if (is.null(dim(data))) {
    
    result <- ApplyFilter(data, filter = (rep(1 / window, window)),
                          method = method, na.rm = TRUE)
    
  } else {
    
    result <- apply(data, 2, ApplyFilter, filter = (rep(1 / window, window)),
                    method = method, na.rm = TRUE)
    
  }
  
  if (hasAgeColumn) {
    result <- data.frame(Year, result)
    colnames(result) <- colNames
  }
    
  return(result)
  
}

#' Estimate isotope slopes
#'
#' Estimate the linear trend slopes in an isotope record within a given
#' moving window.
#'
#' @param x a data frame with two columns with the time scale as the first and
#'   the isotope values as the second column.
#' @param window an integer value giving the size of the moving window within
#'   which the linear trend slopes are estimated; must be an odd number.
#' @return a data frame with two columns 'Year' and 'slope' with the time scale
#'   as the first and the slope values as the second column.
#' @author Thomas Münch
#'
estimateSlopes <- function(x, window = 21) {

  if (window < 3) stop("No sensible estimation window size.")
  if ((window %% 2) == 0) stop("Odd estimation window size expected.")
  
  slopes <- rep(NA, nrow(x))
  window.half <- (window - 1) / 2
  for (i in (window.half + 1) : (nrow(x) - window.half)) {

    j <- (i - window.half) : (i + window.half)

    slopes[i] <- coefficients(lm(x[j, 2] ~ x[j, 1]))[2]

  }

  res <- data.frame(Year = x[, 1], slope = slopes)
  return(res)

}

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

doRunningCorrelation <- function(x, y, window = 101) {

  n <- length(x)

  if ((window %% 2) == 0) stop("Use an odd window size.", call. = FALSE)
  if (window > n) stop("window is larger than length of data.", call. = FALSE)

  if (n != length(y)) stop("x and y must have the same length.", call. = FALSE)

  half <- (window - 1) / 2

  correlation <- rep(NA, n)
  for (i in (half + 1) : (n - half)) {

    subset <- (i - half) : (i + half)
    correlation[i] <- cor(x[subset], y[subset], use = "pair")

  }

  return(correlation)

}
