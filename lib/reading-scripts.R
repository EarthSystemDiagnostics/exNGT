#Function reading data
#Function to read NGT data set
#by Maria Hoerhold and Thomas Muench
#' @title  Read NGT records
#' @description read all NGT and associated d18O records 
#' @param Input path
#' @param read csv.
#' @return dataframe

readNGT<-function(path = "data/NGT2012_AnnualMean_FINAL.csv") {
  if (file.exists(path)) {
    read.csv(path, header = TRUE)
  } else {
    stop("No such file or directory; check path or permissions.")
  }
}


readArctic2k<-function(path = "data/Reconstruction_Arc2kv1.1.1.csv") {
  
  # skip unneeded data columns
  colClasses <- c(rep(NA, 4), rep("NULL", 6))
  
  col.names <- c("Year", "TempAnomaly", "2SigmaLow", "2SigmaHigh")
  
  dat <- read.csv(path, header = TRUE, colClasses = colClasses)
  
  colnames(dat) <- col.names
  
  return(dat)
  
}

processNGT <- function(path = NULL, reference.period = 1990 : 1961) {

  # Read data
  if (is.null(path)) {
    ngt <- readNGT()
  } else {
    ngt <- readNGT(path = path)
  }

  # Skip one record which is not used
  ngt$`B22_12` <- NULL

  # Produce anomaly time series
  for (i in 2 : ncol(ngt)) {
    ngt[, i] <- makeAnomalies(ngt[, i], ngt$Year,
                              reference.period = reference.period)
  }

  return(ngt)

}
