##
## aim:
##   script to read average 60-90 N annual mean HadCrut v5.0.1 netcdf data file
## relation:
##   NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##
## Thomas Muench, AWI, 2021
##

#' Read HadCrut netcdf data and convert to data frame
#'
#' @param path path to your copy of the 'exNGT' repository.
#' @author Thomas Münch
processHadCrut <- function(path) {

  require(ncdf4)
  require(readr)

  infile <- paste0("HadCRUT.5.0.1.0.analysis.anomalies.",
                 "ensemble_mean.annual.60to90_mean.nc")
  infile <- file.path(path, "data-raw", "in-other", infile)

  nc <- ncdf4::nc_open(infile, readunlim = FALSE, suppress_dimvals = TRUE)
  time <- ncdf4::ncvar_get(nc, varid = "time")
  dat <- c(ncdf4::ncvar_get(nc, varid = "tas_mean"))
  ncdf4::nc_close(nc)

  hadcrutArctic <- data.frame(
    Year = seq(1850, length.out = length(time)),
    TempAnomaly = dat)

  outfile <- file.path(path, "data", "HadCrut_5.0.1_Arctic_annual.csv")

  readr::write_csv(hadcrutArctic, file = outfile, na = "")

}
