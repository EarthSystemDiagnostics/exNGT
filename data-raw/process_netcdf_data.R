##
## aim:
##   script to read netcdf data files and process the data into format used for
##   paper analyses
## relation:
##   NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##
## Thomas Muench, AWI, 2021
##

#' Read HadCrut v5.0.1 netcdf data of average 60-90 N annual mean surface
#' temperature time series and convert to data frame.
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

#' Read 20CRv3 netcdf reanalysis data of 50-90 N field of near surface
#' temperature time series, convert to 2D array and save as R data file.
#'
#' @param path path to your copy of the 'exNGT' repository.
#' @author Thomas Münch
process20CR <- function(path) {

  infile <- "noaa.cires.doe.20crv3.air.2m.annual.mean.50to90.nc"
  infile <- file.path(path, "data-raw", "in-other", infile)

  if (!file.exists(infile)) {
    stop(sprintf("%s: Invalid file name.", infile), call. = FALSE)
  }

  nc <- ncdf4::nc_open(infile, readunlim = FALSE, suppress_dimvals = TRUE)

  lon  <- ncdf4::ncvar_get(nc, varid = "lon")
  lat  <- ncdf4::ncvar_get(nc, varid = "lat")
  time <- ncdf4::ncvar_get(nc, varid = "time")
  dat  <- ncdf4::ncvar_get(nc, varid = "air")
  ncdf4::nc_close(nc)

  # time and coordinate information at each grid cell
  time <- seq(1836, length.out = length(time))
  lons <- rep(lon, length(lat))
  lats <- rep(lat, each = length(lon))

  # shape in 2D array
  dim(dat) <- c(length(lat) * length(lon), length(time))
  dat <- t(dat)

  TwenCR <- list(dat = dat, time = time, lat = lats, lon = lons)

  outfile <- file.path(path, "data",
                       "NOAA_CIRES_DOE_20CR_v3_50to90_annual_field.rda")
  save(TwenCR, file = outfile)

}
