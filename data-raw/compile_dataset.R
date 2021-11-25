##
## aim:
##   compile all new and published isotope data and merge it to single file.
## relation:
##   NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##
## Thomas Muench, AWI, 2021
##

#' Compile exNGT 2012 annual isotope data set
#'
#' This function calculates the annual mean isotope data of the new extended NGT
#' 2012 cores (namely, B18_2012, B21_2012, B23_2012, B26_2012, and NGRIP_2012),
#' reads all other published North Greenland annual mean isotope time series,
#' and compiles all data into a single data set which is saved to disk in the
#' ./data folder.
#'
#' @param path path to your copy of the 'exNGT' repository containing the
#'   relevant data and source files needed.
#' @param filename output file name for the compiled data set (any file
#'   extension is ignored); this file will be written in csv format to the
#'   ./data/ folder.
#' @author Thomas Münch
compileDataset <- function(path, filename) {

  require(pangaear)
  require(pdftools)
  require(readr)
  require(stringr)
  require(tibble)
  require(purrr)
  require(plyr)
  require(dplyr)

  source(file.path(path, "data-raw/calculate_year_means_exNGT_2012.R"))
  source(file.path(path, "data-raw/read_pub_data.R"))

  # ----------------------------------------------------------------------------
  # process annual data for the re-drilled NGT cores

  exNGT <- calculateAnnualMeansExNGT2012(path = path)

  # ----------------------------------------------------------------------------
  # read annual data of the original NGT cores + GRIP from PANGAEA repository

  NGT <- getPangaeaData()

  # ----------------------------------------------------------------------------
  # read published annual data of GISP2 and NGRIP cores

  gisp2 <- getGISP2()
  ngrip <- getNGRIP()

  # ----------------------------------------------------------------------------
  # read annual NEGIS and NEEM core data provided by Bo Vinther

  negis <- readr::read_csv(file.path(path, "data-raw/in-other",
                                     "NEGIS_AnnualMean.csv"),
                           col_types = cols())
  neem <- readr::read_csv(file.path(path, "data-raw/in-other",
                                    "NEEM_AnnualMean.csv"),
                          col_types = cols())

  # ----------------------------------------------------------------------------
  # compile to single data frame and save

  filename <- paste0(tools::file_path_sans_ext(filename), ".csv")

  c(exNGT, NGT, list(gisp2 = gisp2, ngrip = ngrip,
                     negis = negis, neem = neem)) %>%
    purrr::reduce(dplyr::full_join, by = "Year") %>%
    readr::write_csv(
      file = file.path(path, "data", filename),
      na = "")
  
}
