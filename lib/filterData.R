#' Filter time series
#'
#' Apply a given filter to a time series using different endpoint constraints as
#' described below.
#'
#' Note that when passing objects of class \code{ts}, the time step provided is
#' not used; thus, for time series with a time step different from 1, the filter
#' has to be adapted accordingly.
#'
#' The function applies endpoint constrains following Mann et al., GRL, 2003.
#' Available methods are:
#' method = 0: no constraint (loss at both ends);
#' method = 1: minimum norm constraint;
#' method = 2: minimum slope constraint;
#' method = 3: minimum roughness constraint;
#' method = 4: circular filtering.
#'
#' @param data numeric vector with the input timeseries (ts object).
#' @param filter numeric vector of filter weights.
#' @param method single integer for choosing an endpoint constraint method;
#'   available choices are integers 0-4.
#' @return filtered timeseries (ts object).
#' @author Thomas Laepple
#'
ApplyFilter <- function(data, filter, method = 0) {

  if (!method %in% (0 : 4))
    stop("Unknown method; only 0 : 4 available.")

  circular = FALSE

  if (method == 0 | method == 4) {

    if (method == 4) {circular = TRUE}

    result <- stats::filter(c(data), filter, circular = circular)

  } else {

    N <- floor(length(filter) / 2)

    if (method == 1) {

      before <- rep(mean(data), N)
      after  <- rep(mean(data), N)

    } else if (method == 2) {

      before <- c(data)[N : 1]
      after  <- c(data)[length(data) : (length(data) - N + 1)]

    } else if (method == 3) {

      before <- c(data)[N : 1]
      after  <- c(data)[length(data) : (length(data) - N + 1)]
      
      before <- c(data)[1] - (before - mean(before))
      after  <- c(data)[length(data)] - (after - mean(after))

    }

    result <- stats::filter(c(before, data, after), filter, circular = circular)
    result <- result[(N + 1) : (N + length(data))]
  }

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
  
  if (hasAgeColumn) {
    colNames <- colnames(data)
    Year <- data[, 1]
    data <- data[, -1]
  }
  
  if (is.null(dim(data))) {
    
    result <- ApplyFilter(data, filter = (rep(1 / window, window)),
                          method = method)
    
  } else {
    
    result <- apply(data, 2, ApplyFilter, filter = (rep(1 / window, window)),
                    method = method)
    
  }
  
  if (hasAgeColumn) {
    result <- data.frame(Year, result)
    colnames(result) <- colNames
  }
    
  return(result)
  
}
