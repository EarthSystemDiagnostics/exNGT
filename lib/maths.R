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
#'
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

#' Calculate running correlation
#'
#' Calculate the running correlation between two time series over consecutive
#' windows of a given length.
#'
#' @param x numeric vector (time series) to be correlated with \code{y}.
#' @param y numeric vector (time series) to be correlated with \code{x}.
#' @param window integer; width of the windows in number of time steps or
#'   indices over which correlations between \code{x} and \code{y} are computed;
#'   needs to be an odd number.
#' @return numeric vector of the same length as \code{x} and \code{y} with the
#'   running correlations between the data computed over every consecutive
#'   window which fits into the observational range of the data. Note that the
#'   correlation estimates are centred on the window, so that the first and last
#'   \code{(window - 1) / 2} return values are \code{NA}.
#' @author Thomas Münch
#'
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

#' Estimate correlation and its significance
#'
#' Calculate the Pearson correlation coefficient between two running-mean
#' filtered time series and estimate its significance based on a Monte Carlo
#' sampling of surrogate data.
#'
#' The significance of the correlation between the data and the reference signal
#' is estimated as follows: The AR1 autocorrelation coefficient is estimated
#' from the original (i.e. unfiltered) data and used to create random surrogate
#' time series with the same autocorrelation structure. The surrogate time
#' series are filtered with the same filter as the data, and for each of the
#' filtered surrogate time series the correlation with the filtered signal is
#' computed. The p value of the correlation between the data and the signal is
#' then obtained from the fraction of surrogate correlation values which lie
#' above the actual observed correlation between signal and data.
#'
#' @param data a data frame with the unfiltered data time series (with the time
#'   points in the first column and the data in the second column). This is only
#'   used to estimate the autocorrelation in the data, so in order to get a
#'   reliable autocorrelation estimate it may cover a larger time interval than
#'   the actual \code{filteredData} (assuming stationarity).
#' @param filteredData a data frame with the filtered version of \code{data},
#'   which is correlated with the filtered \code{signal}.
#' @param signal a data frame with the filtered signal time series (with the
#'   time points in the first column and the data in the second column).
#' @param filter.window single integer giving the size of the running mean
#'   filter window for filtering the surrogate data; needs to be identical to
#'   the one applied for filtering the input data.
#' @param analysis.period optional vector of time points to define the time
#'   period over which the correlation is estimated; if \code{NULL} (the
#'   default) all time points which have non-missing observations in both the
#'   \code{filteredData} and the \code{signal} are used.
#' @param timeColumn character vector of length 1 with the name of the time
#'   column for \code{filteredData} and \code{signal}.
#' @param nmc integer; number of surrogate time series to use for estimating the
#'   correlation significance.
#' @return a list of two elements: the correlation between the
#'   \code{filteredData} and the \code{signal} and the estimated p value.
#' @author Thomas Münch
#'
estimateCorrelation <- function(data, filteredData, signal,
                                filter.window, analysis.period = NULL,
                                timeColumn = "Year", nmc = 1000) {

  a1 <- acf(data[, 2], lag.max = 1, plot = FALSE)$acf[2, 1, 1]

  if (length(analysis.period)) {
    filteredData <- filteredData %>%
      dplyr::filter(!!as.name(timeColumn) %in% analysis.period)
  }

  corData <- filteredData %>%
    dplyr::inner_join(signal, by = timeColumn) %>%
    na.omit()

  referenceCor <- cor(corData[, 2], corData[, 3])

  n <- nrow(corData)
  surrogates <- replicate(nmc, arima.sim(model = list(ar = a1), n = n))

  filteredSurrogates <- filterData(surrogates, window = filter.window,
                                   hasAgeColumn = FALSE)

  surrogateCor <- cor(corData[, 2], filteredSurrogates)[1, ]

  p <- sum(surrogateCor >= referenceCor) / nmc

  return(list(r = referenceCor, p = p))

}

#' Estimate running correlation and its significance
#'
#' Calculate the running correlation between two time series over consecutive
#' windows of a given length and estimate the significance of the correlations
#' based on a Monte Carlo sampling of surrogate data, accounting for running
#' mean filtering of the data.
#'
#' The significance of the running correlation between the data and the
#' reference signal is estimated following the method described in van
#' Oldenborgh and Burgers (2005): The overall correlation between the data and
#' the signal is computed and used to create random surrogate time series which
#' on average exhibit the same correlation with the signal. The data, signal and
#' surrogate data are then filtered with a running mean filter, and the running
#' correlation between the data and the signal and between the surrogate data
#' and the signal are computed. Expressing the correlations in terms of z values
#' (van Oldenborgh and Burgers, 2005), the p value of the running correlation is
#' obtained from the fraction of maximum z value differences for the surrogate
#' data which exceed the maximum z value difference of the observation.
#'
#' @param data numeric vector with the investigated time series that is to be
#'   correlated with the reference \code{signal}; at the original temporal
#'   resolution prior to any running mean filtering.
#' @param signal numeric vector of the reference signal time series to which the
#'   \code{data} is correlated; at the same temporal resolution as \code{data}.
#' @param nmc integer; number of surrogate time series to use for estimating the
#'   correlation significance.
#' @param correlation.window integer; width of the windows in number of time
#'   steps or indices over which correlations between \code{data} and
#'   \code{signal} are computed; needs to be an odd number (see also
#'   \code{doRunningCorrelation()}).
#' @param filter.window single integer giving the size of the running mean
#'   filter window for filtering the data.
#' @return a list of six elements: estimated p value, time series of running
#'   correlation between signal and data, array of running correlation time
#'   series between signal and surrogate data, arrays with the time series of
#'   the 2.5 and 97.5 % quantiles of the running correlations between signal and
#'   surrogate data, and execution date of the analysis.
#' @references
#' van Oldenborgh and Burgers, Searching for decadal variations in ENSO
#'   precipitation teleconnections, Geophys. Res. Lett., 32(15), L15701,
#'   https://doi.org/10.1029/2005GL023110, 2005.
#' @author Thomas Münch
#'
estimateRunningCorrelation <- function(data, signal, nmc = 1000,
                                       correlation.window = 101,
                                       filter.window = 11) {

  deltaZ <- function(r) {

    z <- function(r) {log((1 + r) / (1 - r)) / 2}

    z(max(r, na.rm = TRUE)) - z(min(r, na.rm = TRUE))
  }

  r    <- cor(data, signal)
  rfac <- sqrt((1 - r^2) / r^2)

  surrogates <- replicate(nmc, c(scale(signal)) + rfac * rnorm(length(signal)))

  signal <- filterData(signal, window = filter.window, hasAgeColumn = FALSE)
  data <- filterData(data, window = filter.window, hasAgeColumn = FALSE)
  surrogates <- apply(surrogates, 2, filterData, window = filter.window,
                      hasAgeColumn = FALSE)

  runningDat <- doRunningCorrelation(data, signal, correlation.window)
  runningSim <- apply(surrogates, 2, doRunningCorrelation, y = signal,
                      window = correlation.window)

  deltaZDat <- deltaZ(runningDat)
  deltaZSim <- apply(runningSim, 2, deltaZ)

  res <- list(
    p = sum(deltaZSim >= deltaZDat) / nmc,
    dat = runningDat,
    sim = runningSim,
    r.upper = apply(runningSim, 1, quantile, probs = 0.975, na.rm = TRUE),
    r.lower = apply(runningSim, 1, quantile, probs = 0.025, na.rm = TRUE),
    date = Sys.time()
  )

  return(res)

}

#' Spectral coherence estimate
#'
#' This functions calculates the magnitude-squared coherence between two time
#' series using the smoothed periodogram with confidence level estimated from
#' a Monte Carlo procedure.
#'
#' @param x1 numerical vector with timeseries one.
#' @param x2 numerical vector with timeseries two of the same length as
#'   \code{x1}.
#' @param spans vector of odd integers giving the widths of modified Daniell
#'   smoothers to be used to smooth the periodogram.
#' @param nmc number of Monte Carlo realizations used for estimating the
#'   confidence level.
#' @param p significance threshold (p value) used for the confidence level.
#' @return a list of three elements: frequency vector and vector of associated
#'   magnitude-squared coherence estimates and the global confidence level
#'   (single number).
#' @author Thomas Laepple
#'
coherence <- function(x1, x2, spans, nmc = 100, p = 0.95) {

  if (length(x1) != length(x2)) {
    stop("'x1' and 'x2' need to have the same length.")
  }

  coh   <- spectrum(cbind(x1, x2), spans = spans, plot = FALSE)
  cohCL <- mcCoherence(x1, x2, spans = spans, nmc = nmc, p = p)

  return(list(freq = coh$freq, coh = c(coh$coh), confLevel = cohCL))

}

#' Confidence level for spectral coherence
#'
#' This functions estimates the confidence level for the magnitude-squared
#' coherence between two time series by replacing the time series with red noise
#' surrogates, calculating the coherence between these surrogates and returning
#' the desired quantile of the resulting coherence estimates as the confidence
#' level.
#'
#' @param x1 numerical vector with timeseries one.
#' @param x2 numerical vector with timeseries two of the same length as
#'   \code{x1}.
#' @param spans vector of odd integers giving the widths of modified Daniell
#'   smoothers to be used to smooth the periodogram.
#' @param nmc number of Monte Carlo realizations (i.e. the number of created
#'   red-noise surrogate pairs) used for estimating the confidence level.
#' @param p significance threshold (p value) used for the confidence level.
#' @return the confidence level for the squared coherence based on the desired
#'   significance value.
#' @author Thomas Laepple
#'
mcCoherence <- function(x1, x2, spans, nmc = 100, p = 0.95) {

  if (length(x1) != length(x2)) {
    stop("'x1' and 'x2' need to have the same length.")
  }

  getA1 <- function(x) {
    acf(c(x), plot = FALSE)$acf[2]
  }

  red <- function(a1, n) {
    c(arima.sim(list(ar = a1), n))
  }

  x1.a1 <- getA1(x1)
  x2.a1 <- getA1(x2)

  foo <- spectrum(x1, plot = FALSE)
  mcCohMatrix <- matrix(NA, nmc, length(foo$freq))

  mcCohMatrix <- sapply(1 : nmc, function(n) {
    spectrum(cbind(red(x1.a1, length(x1)), red(x2.a1, length(x2))),
             spans = spans, plot = FALSE)$coh
  })

  # return mean confidence level since level is independent of frequency
  mean(apply(mcCohMatrix, 2, quantile, probs = p))

}
