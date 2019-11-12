##' Anomaly time series
##'
##' Calculate anomaly time series relative to a given reference period.
##' @param data numeric vector (time series)
##' @param age the age scale of \code{data}
##' @param reference.period numeric vector of ages of the reference period
##' @return numeric vector with the anomaly time series
##' @author Thomas Münch
makeAnomalies <- function(data, age, reference.period = 1990 : 1961) {

  data - mean(data[match(reference.period, age)], na.rm = TRUE)
}
