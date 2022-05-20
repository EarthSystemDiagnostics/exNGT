##
## aim:
##   library functions to compile available NGT accumulation data
## relation:
##   NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##
## Thomas Muench, AWI, 2022
##

#' Wrapper function to read in accumulation record
#'
#' @param path path relative to the 'exNGT' repository of the data folder
#'   containing the accumulation record.
#' @param file file name of the accumulation record.
#' @param name firn core name of the accumulation record.
#' @param delim column delimiter used in the data file.
#' @return data frame of two columns: Year CE and accumulation data.
#' @author Thomas Münch
readAccumulationRecord <- function(path, file, name = "", delim = "\t") {

  readr::read_delim(file.path(path, file), col_types = readr::cols(),
                    delim = delim) %>%
    setNames(c("Year", name))

}

#' Read in NEEM accumulation record
#'
#' This function reads in the four available NEEM accumulation records of the
#' firn cores NEEM2010S2, NEEM2008S3, NEEM2007S3 and NEEM2008S2, and calculates
#' the average accumulation record.
#'
#' @param path path relative to the 'exNGT' repository of the data folder
#'   containing the NEEM accumulation records.
#' @param file file name of the excel file with the NEEM accumulation records.
#' @param name record name to use for the averaged record, defaults to "NEEM".
#' @return data frame of two columns: Year CE and average accumulation data.
#' @author Thomas Münch
readNeemAccumulation <- function(path = "data-raw/in-other",
                                 file = "NEEM_Acc.xlsx", name = "NEEM") {

  tmp <- list(readxl::read_xlsx(path = file.path(path, file),
                                range = "A2:B157") %>%
              setNames(c("Year", "NEEM2010S2")),
              readxl::read_xlsx(path = file.path(path, file),
                                range = "D2:E284") %>%
              setNames(c("Year", "NEEM2008S3")),
              readxl::read_xlsx(path = file.path(path, file),
                                range = "G2:H267") %>%
              setNames(c("Year", "NEEM2007S3")),
              readxl::read_xlsx(path = file.path(path, file),
                                range = "J2:K156") %>%
              setNames(c("Year", "NEEM2008S2"))) %>%
    purrr::reduce(dplyr::full_join, by = "Year") %>%
    as.data.frame()

  data.frame(Year = tmp$Year, NEEM = rowMeans(tmp[, -1], na.rm = TRUE))

}

#' Compile NGT 2012 accumulation data set
#'
#' This function reads in all available NGT accumulation data, compiles them
#' into a single data frame, and saves the data as a csv file.
#'
#' @param path path to your copy of the 'exNGT' repository.
#' @param filename output file name for the compiled data set (any file
#'   extension is ignored); this file will be written in csv format to the
#'   ./data/ folder.
#' @author Thomas Münch
compileAccumulationRecords <- function(path, filename) {

  require(readr)
  require(readxl)
  require(purrr)
  require(dplyr)
  
  path1 <- "data-raw/out-pangaea"
  path2 <- "data-raw/in-other"

  filename <- paste0(tools::file_path_sans_ext(filename), ".csv")

  # read all accumulation records, compile into single data frame, and save
  list(
    b18.12 = readAccumulationRecord(path1, file = "B18_2012_AccmRate.txt",
                                    name = "B18_12"),
    b21.12 = readAccumulationRecord(path1, file = "B21_2012_AccmRate.txt",
                                    name = "B21_12"),
    b23.12 = readAccumulationRecord(path1, file = "B23_2012_AccmRate.txt",
                                    name = "B23_12"),
    b26.12 = readAccumulationRecord(path1, file = "B26_2012_AccmRate.txt",
                                    name = "B26_12"),
    ngrip.12 = readAccumulationRecord(path1, file = "NGRIP_2012_AccmRate.txt",
                                      name = "NGRIP_12", delim = " "),
    b16 = readAccumulationRecord(path2, file = "B16_Schwager_Accmrate_we.txt",
                                 name = "B16"),
    b18 = readAccumulationRecord(path2, file = "B18_Schwager_Accmrate_we.txt",
                                 name = "B18"),
    b21 = readAccumulationRecord(path2, file = "B21_Schwager_Accmrate_we.txt",
                                 name = "B21"),
    b26 = readAccumulationRecord(path2, file = "B26_Schwager_Accmrate_we.txt",
                                 name = "B26"),
    b29 = readAccumulationRecord(path2, file = "B29_Schwager_Accmrate_we.txt",
                                 name = "B29"),
    neem = readNeemAccumulation()
       ) %>%
    purrr::reduce(dplyr::full_join, by = "Year") %>%
    dplyr::filter(Year <= 2011) %>%
    dplyr::mutate(dplyr::across(!Year, ~ .x * 1.e3)) %>% # in mm w.eq.
    as.data.frame() %>%
    round() %>% # round to mm resolution
    readr::write_csv(file = file.path(path, "data", filename), na = "")

}
