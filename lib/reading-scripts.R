#' Read NGT records
#'
#' Read all NGT and associated oxygen isotope records.
#'
#' @param path file path (relative to working directory) of the NGT data set.
#' @return a data frame with the read NGT data set.
#' @author Thomas Münch, Maria Hoerhold
#'
readNGT <- function(path = "data/NGT2012_AnnualMean_FINAL.csv") {

  if (file.exists(path)) {
    read.csv(path, header = TRUE)
  } else {
    stop("No such file or directory; check path or permissions.")
  }
}

#' Read Arctic2k temperature reconstruction
#'
#' Read the annual Arctic2k temperature reconstruction record of temperature
#' anomalies w.r.t. 1961-1990 including +/- 2 sigma uncertainty.
#'
#' @param path file path (relative to working directory) of the Arctic2k data
#'   set.
#' @return a data frame of four columns and 2000 rows with the read Arctic2k
#'   data set.
#' @author Thomas Münch, Maria Hoerhold
#'
readArctic2k <- function(path = "data/Reconstruction_Arc2kv1.1.1.csv") {
  
  # skip unneeded data columns
  colClasses <- c(rep(NA, 4), rep("NULL", 6))
  
  colNames <- c("Year", "TempAnomaly", "2SigmaLow", "2SigmaHigh")
  
  dat <- read.csv(path, header = TRUE, colClasses = colClasses)
  
  colnames(dat) <- colNames
  
  return(dat)
  
}

#' Process NGT records
#'
#' Process the NGT and associated oxygen isotope records, which includes reading
#' in the data, skipping unsuitabe records, and calculating anomaly time series.
#'
#' @param path file path of the NGT data set. For \code{NULL} (the default) the
#'   default file path set by the reading function \code{readNGT} is used.
#' @param reference.period numeric vector of ages of the reference period for
#'   calculating the anomaly time series; defaults to the standard 1961-1990
#'   period.
#' @return a data frame of anomaly time series of the relevant NGT and
#'   associated isotope records.
#' @author Thomas Münch
#'
processNGT <- function(path = NULL, reference.period = 1990 : 1961) {

  # Read data
  if (is.null(path)) {
    ngt <- readNGT()
  } else {
    ngt <- readNGT(path = path)
  }

  # Skip two records not used
  ngt$`B22_12` <- NULL # only has data from 1997-2011
  ngt$`B19`    <- NULL # has no data in reference period 1961-1990

  # Produce anomaly time series
  ngt <- makeAnomalies(ngt)

  return(ngt)

}

#' Read DMI instrumental temperatures
#'
#' Read the Danish Meteorological Insitute (DMI) instrumental temperature
#' records of the Greenlandic stations Pituffik, Upernavik and Danmarkshavn.
#'
#' @param path file path (relative to working directory) of the DMI data
#'   set.
#' @return a data frame of four columns and 146 rows with the read DMI
#'   data set for the three stations.
#' @author Thomas Münch
#'
readDMI <- function(path = "data/gr_annual_temperature_1873_2015.csv") {

  colNames <- c("Year", "Pituffik", "Upernavik", "Danmarkshavn")

  dat <- read.csv(path, header = TRUE, na.strings = "#NULL!")
  dat <- dat[nrow(dat) : 1, c(1, 2, 4, 12)]

  colnames(dat) <- colNames

  return(dat)

}
