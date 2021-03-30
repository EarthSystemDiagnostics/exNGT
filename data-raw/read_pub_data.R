##
## aim:
##   library functions to read published part of the used Greenland isotope data
## relation:
##   NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##
## Thomas Muench, AWI, 2021
##

#' Download and parse NGT and GRIP annual isotope data from PANGAE repository
#'
#' @return a list of tibbles for each isotope data set with the first column
#'   being time in years CE and the second column being the annual isotope
#'   value in permil.
#' @author Thomas Münch
getPangaeaData <- function() {

  require(pangaear)
  require(tibble)
  require(dplyr)
  
  pangaeaInfo <- tibble::tribble(
    ~Line, ~`Site`,  ~`DOI`,                   ~`dataCol`,
    # -- / ------ /  ----- /                   --------- /
    1,     "B16",    "10.1594/PANGAEA.849148", 6,
    2,     "B17",    "10.1594/PANGAEA.849149", 6,
    3,     "B18",    "10.1594/PANGAEA.849150", 6,
    4,     "B20",    "10.1594/PANGAEA.849152", 6,
    5,     "B21",    "10.1594/PANGAEA.849153", 6,
    6,     "B22",    "10.1594/PANGAEA.849154", 6,
    7,     "B23",    "10.1594/PANGAEA.849155", 6,
    8,     "B26",    "10.1594/PANGAEA.849156", 6,
    9,     "B27.28", "10.1594/PANGAEA.849160", 3,
    10,    "B29",    "10.1594/PANGAEA.849237", 6,
    11,    "B30",    "10.1594/PANGAEA.849159", 6,
    12,    "GRIP",   "10.1594/PANGAEA.786354", 5
  )

  pangaeaData <- list()
  for (i in 1 : nrow(pangaeaInfo)) {

    pangaeaData[[i]] <- pangaeaInfo$DOI[i] %>%
      pangaear::pg_data() %>%
      .[[1]] %>%
      .$data %>%
      dplyr::select("Age [a AD]", pangaeaInfo$dataCol[i]) %>%
      setNames(c("Year", pangaeaInfo$Site[i]))
  }

  return(setNames(pangaeaData, pangaeaInfo$Site))

}

#' Download and parse GISP2 annual isotope data from U Washington archive
#'
##' @return a tibble with two variables: time in years CE and annual isotope
#'   value in permil.
#' @author Thomas Münch
getGISP2 <- function() {

  require(readr)
  require(dplyr)
  
  file <- "https://depts.washington.edu/qil/datasets/gisp2_1yr.txt"

  readr::read_delim(file, delim = " ", skip = 12,
                    col_names = c("depth", "Year", "GISP2"),
                    col_types = "nnn") %>%
    dplyr::mutate(Year = 1950 - Year) %>%
    dplyr::select(-depth) %>%
    round(2)

}

#' Download and parse NGRIP annual isotope data from U Copenhagen archive
#'
#' @param tmpdir directory in which to temporarily store the downloaded data
#'   file; the file will be deleted from there once the function has parsed the
#'   data.
#' @return a tibble with two variables: time in years CE and annual isotope
#'   value in permil.
#' @author Thomas Münch
getNGRIP <- function(tmpdir = "/tmp") {

  require(pdftools)
  require(readr)
  require(stringr)
  require(plyr)
  require(dplyr)

  con <- "https://www.iceandclimate.nbi.ku.dk/data"
  file <- "Kaufman_etal_2009_data_29sep2009.pdf"

  file_path <- file.path(tmpdir, file)

  download.file(url = file.path(con, file), destfile = file_path)

  raw <- pdftools::pdf_text(pdf = file_path) %>%
    readr::read_lines() %>%    # read line by line
    .[-(1 : 31)] %>%           # remove intro lines
    stringr::str_squish() %>%  # remove repeated whitespace
    .[-which(nchar(.) == 0)]   # remove zero length lines

  # there are peculiar pdf lines splitted over sets of four lines
  i <- which(nchar(raw) < 10)
  j <- i[c(which(diff(i) > 1) - 3, length(i) - 3)]

  # concatenate the proper with the peculiar lines
  cleaned <- character()
  for (index in 1 : length(j)) {

    i.start <- 1
    if (index != 1) i.start <- j[index - 1] + 4

    start <- raw[i.start : (j[index] - 1)]
    end   <- paste0(raw[j[index] : (j[index] + 3)], collapse = "")

    # remove the extra commas
    x <- which(strsplit(end, "")[[1]] == ",")[c(1, 3, 5)]
    offset <- 0
    for (ix in 1 : length(x)) {
      stringr::str_sub(end, x[ix] - offset, x[ix] - offset) <- ""
      offset <- offset + 1
    }

    cleaned <- c(cleaned, start, end)
  }
  cleaned <- c(cleaned, raw[(j[index] + 4) : length(raw)])

  # replace decimal character, convert time axis,
  # and combine relevant data into tibble
  parsed <- cleaned %>%
    stringr::str_replace_all(",", ".") %>%
    strsplit(split = " ") %>%
    plyr::ldply() %>%
    setNames(c("Year", "NGRIP", "DYE3", "Agassiz")) %>%
    transform(Year = as.numeric(Year), NGRIP = as.numeric(NGRIP)) %>%
    dplyr::mutate(Year = 2000 - Year) %>%
    dplyr::select(Year, NGRIP) %>%
    dplyr::as_tibble()

  # remove input file
  system(sprintf("rm %s", file_path))

  return(parsed)

}
