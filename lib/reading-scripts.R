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
