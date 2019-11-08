#Function reading data
#Function to read NGT data set
#by Maria Hoerhold and Thomas Muench
#' @title  Read NGT records
#' @description read all NGT and associated d18O records 
#' @param Input path
#' @param read csv.
#' @return dataframe

readNGT<-function(path = "data/NGT2012_AnnualMean_FINAL.csv") {
  read.csv(path, header = TRUE)
}


readArctic2k<-function(path = "data/Reconstruction_Arc2kv1.1.1.csv") {
  
  # skip unneeded data columns
  colClasses <- c(rep(NA, 4), rep("NULL", 6))
  
  col.names <- c("Year", "TempAnomaly", "2SigmaLow", "2SigmaHigh")
  
  dat <- read.csv(path, header = TRUE, colClasses = colClasses)
  
  colnames(dat) <- col.names
  
  return(dat)
  
}