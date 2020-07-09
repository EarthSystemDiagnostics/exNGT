#' Anomaly time series
#'
#' Calculate anomaly time series relative to the mean of a given reference
#' period.
#'
#' @param data a data set for which to calculate the anomalies; can be a single
#'   numeric vector or a data frame or matrix of vectors. In the latter case,
#'   the first column can be used to provide the age scale of the data; else use
#'   the \code{age} parameter for this.
#' @param age provide the age scale of \code{data} here if it is not supplied
#'   in the latter parameter itself.
#' @param reference.period numeric vector with the ages of the reference period.
#' @param hasAgeColumn logical; do the supplied \code{data} contain an age
#'   column as the first column? Setting this to \code{FALSE} but leaving
#'   \code{age} as \code{NULL} will cause an error.
#' @return numeric vector, matrix or data frame with the anomaly time series.
#' @author Thomas Münch
#'
makeAnomalies <- function(data, age = NULL, reference.period = 1990 : 1961,
                          hasAgeColumn = TRUE) {

  if (is.null(age) & !hasAgeColumn) {
    stop("Age scale information missing.", call. = FALSE)
  }

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
  } else {
    Year <- age
  }

  i <- match(reference.period, Year)

  if (is.null(dim(data))) {

    result <- data - mean(data[i], na.rm = TRUE)

  } else {

    result <- apply(data, 2, function(x) {
      x - mean(x[i], na.rm = TRUE)})

  }

  if (hasAgeColumn) {
    result <- data.frame(Year, result)
    colnames(result) <- colNames

  }

  return(result)
}
