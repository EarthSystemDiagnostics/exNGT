##
## aim:
##   calculate the annual means of the NGT 2012 extension firn core records
##   based on the depth-age relationships for each core.
## relation:
##   NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##
## Thomas Muench, AWI, 2021
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
  require(PaleoSpec)

  path <- file.path(path, "data-raw")

  # ----------------------------------------------------------------------------
  # read age-depth relationship and raw isotope depth profiles

  sites <- c("B18", "B21", "B23", "NGRIP")

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

    dated[[site]] <- PaleoSpec::AvgToBin(x = raw[[site]]$depth_center,
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

