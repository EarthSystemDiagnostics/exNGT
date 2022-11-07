##
## aim:
##   calculate the annual means of the NGT 2012 extension firn core records
##   based on the depth-age relationships for each core.
## relation:
##   NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##
## Thomas Muench, AWI, 2022
##

#' Calculate exNGT 2012 annual means
#'
#' Calculate the annual mean d18O values for the extended NGT 2012 cores, based
#' on the raw isotope depth profiles and the age-depth relationships.
#'
#' @param path path to your copy of the 'exNGT' repository containing the
#'   relevant data needed.
#' @param update logical; shall the processed annual mean data be written to
#'   disk, which overrides (updates) existing files?
#' @return invisibly the processed annnual mean data as a list.
#' @author Thomas Münch
calculateAnnualMeansExNGT2012 <- function(path, update = FALSE) {

  require(readr)
  require(dplyr)
  require(tibble)
  require(purrr)

  path <- file.path(path, "data-raw")

  # ----------------------------------------------------------------------------
  # read age-depth relationship and raw isotope depth profiles

  sites <- c("B18", "B21", "B23", "NGRIP", "B26")

  files <- sprintf("%s/in-raw/%s_2012_DepthAge.txt", path, sites)

  dat <- files %>%
    lapply(function(f) {
      x <- f %>%
        readr::read_tsv(col_types = readr::cols()) %>%
        setNames(c("Year", "depth_top", "depth_bottom"))}) %>%
    setNames(sites)

  files <- sprintf("%s/in-raw/%s_2012_IsotopeData.txt", path, sites)

  raw <- files %>%
    lapply(function(f) {
      x <- f %>%
        readr::read_tsv(col_types = readr::cols()) %>%
        setNames(c("depth_center", "d18O"))}) %>%
    setNames(sites)

  # ----------------------------------------------------------------------------
  # calculate annual means

  dated <- list()
  for (site in sites) {

    dated[[site]] <- AvgToBin(x = raw[[site]]$depth_center,
                              y = raw[[site]]$d18O,
                              breaks = dat[[site]]$depth_top,
                              right = FALSE) %>%
      .[-1] %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Year = rev(rev(dat[[site]]$Year)[-1])) %>%
      dplyr::select(Year, avg) %>%
      setNames(c("Year", sprintf("%s_12", site))) %>%
      dplyr::filter(`Year` != 2012)

  }

  # ----------------------------------------------------------------------------
  # output files for PANGAEA

  if (update) {

    purrr::walk(names(dated), function(nm) {
      dated[[nm]] %>%
        round(2) %>%
        readr::write_tsv(
          file = sprintf("%s/out-pangaea/%s_2012_AnnualMean.txt", path, nm),
          na = "")
    })
  }

  dated %>%
    lapply(function(x) {round(x, 2)}) %>%
    invisible()

}

#' Calculate NEGIS annual means
#'
#' Calculate the annual mean d18O values for the NEGIS core, based on the
#' published raw isotope depth profile and age-depth relationship.
#'
#' @return invisibly the processed annnual mean data as a data frame.
#' @author Thomas Münch
calculateAnnualMeansNEGIS <- function() {

  require(readr)
  require(dplyr)
  require(tibble)

  # ----------------------------------------------------------------------------
  # read published age-depth relationship and raw isotope depth profile

  con <- "https://www.ncei.noaa.gov/pub/data/paleo/icecore/greenland/"

  dat <- paste0(con, "negis2012age.txt") %>%
    readr::read_tsv(skip = 104, col_types = readr::cols()) %>%
    setNames(c("depth_bottom", "thickness", "Year", "YearBP")) %>%
    dplyr::mutate(depth_top = depth_bottom - thickness) %>%
    dplyr::select(Year, depth_top, depth_bottom) %>%
    dplyr::slice(-1)

  dat <- dat %>%
    dplyr::add_row(Year = rev(dat$Year)[1] - 1,
                   depth_top = rev(dat$depth_bottom)[1],
                   depth_bottom = NA)

  raw <-   paste0(con, "negis2012chem.txt") %>%
    read.table(header = TRUE, skip = 107, sep = "\t") %>%
    setNames(c("depth_center", "d18O", "dxs", "dust", "cond", "na", "nh4")) %>%
    dplyr::select(depth_center, d18O)

  # ----------------------------------------------------------------------------
  # calculate annual means

  dated <- AvgToBin(x = raw$depth_center,
                    y = raw$d18O,
                    breaks = dat$depth_top,
                    right = FALSE) %>%
    .[-1] %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Year = rev(rev(dat$Year)[-1])) %>%
    dplyr::select(Year, avg) %>%
    setNames(c("Year", "NEGIS"))

  dated %>%
    round(2) %>%
    invisible()

}

#' Bin averaging
#'
#' Average a vector into bins.
#'
#' This function averages the vector \code{y} into bins according to the positon
#' of \code{x} within the breaks. You can either specify a desired number N of
#' breaks which are used to calculate the actual breaks via \code{pretty(x, N)},
#' or directly specify the N + 1 break positions. For \code{right = TRUE} (the
#' default) the averaging bins are defined via \code{x > breaks[i]} and \code{x
#' <= breaks[i + 1]}, else they are defined via \code{x >= breaks[i]} and
#' \code{x < breaks[i + 1]}. If \code{bFill = TRUE}, empty bins are filled using
#' linear interpolation from the neighbours to the center of the bin.
#'
#' Probably the binning could be considerably speeded up by using \code{?cut}.
#'
#' @param x vector of values on which the data in \code{y} is tabulated;
#'   e.g. depth or time points.
#' @param y vector of observation values to be averaged into bins. Must have the
#'   same length as \code{x}.
#' @param N desired number of breaks (ignored if \code{breaks} are supplied
#'   directly).
#' @param breaks vector of break point positions to define the averagig bins; if
#'   omitted, break point positions are calculated from the range of \code{x}
#'   and the desired number of breaks given by \code{N}.
#' @param right logical; indicate whether the bin intervals should be closed on
#'   the right and open on the left (\code{TRUE}, the default), or vice versa
#'   (\code{FALSE}).
#' @param bFill logical; if \code{TRUE}, fill empty bins using linear
#'   interpolation from the neighbours to the center of the bin.
#'
#' @return a list with four elements:
#' \describe{
#' \item{\code{breaks}:}{numeric vector of the used break point positions.}
#' \item{\code{centers}:}{numeric vector with the positions of the bin centers.}
#' \item{\code{avg}:}{numeric vector with the bin-averaged values.}
#' \item{\code{nobs}:}{numeric vector with the number of observations
#'   contributing to each bin average.}
#' }
#'
#' @author Thomas Laepple
#' @source https://github.com/EarthSystemDiagnostics/paleospec
AvgToBin <- function(x, y, N = 2, breaks = pretty(x, N),
                     right = TRUE, bFill = FALSE) {

  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }

  nBins <- length(breaks) - 1

  centers <- (breaks[1 : nBins] + breaks[2 : (nBins + 1)]) / 2
  nObs <- avg <- rep(NA, nBins)

  for (i in 1 : nBins) {

    if (right) {
      selection <- y[which((x > breaks[i]) & (x <= breaks[i + 1]))]
    } else {
      selection <- y[which((x >= breaks[i]) & (x < breaks[i + 1]))]
    }

    avg[i]  <- mean(na.omit(selection))
    nObs[i] <- sum(!is.na(selection))

  }

  if ((sum(missing <- is.na(avg)) > 0) & (bFill)) {

    avg[missing] <- (approx(x, y, centers)$y)[missing]

  }

  list(breaks = breaks, centers = centers, avg = avg, nobs = nObs)

}

