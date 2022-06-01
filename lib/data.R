#' Read NGT records
#'
#' Read all NGT and associated oxygen isotope records.
#'
#' @param path file path (relative to working directory) of the NGT data set.
#' @return a data frame with the read NGT data set.
#' @author Thomas Münch, Maria Hoerhold
#'
readNGT <- function(path = "data/NGT2012_AnnualMean_20211124.csv") {

  if (file.exists(path)) {
    read.csv(path, header = TRUE)
  } else {
    stop("No such file or directory; check path or permissions.")
  }
}

#' Read Arctic2k temperature reconstruction
#'
#' Read the annual Arctic2k temperature reconstruction record of temperature
#' anomalies w.r.t. 1961-1990 including +/- 2 sigma uncertainty.
#'
#' @param path file path (relative to working directory) of the Arctic2k data
#'   set.
#' @return a data frame of four columns and 2000 rows with the read Arctic2k
#'   data set.
#' @author Thomas Münch, Maria Hoerhold
#'
readArctic2k <- function(path = "data/Reconstruction_Arc2kv1.1.1.csv") {
  
  # skip unneeded data columns
  colClasses <- c(rep(NA, 4), rep("NULL", 6))
  
  read.csv(path, header = TRUE, colClasses = colClasses) %>%
    setNames(c("Year", "TempAnomaly", "2SigmaLow", "2SigmaHigh")) %>%
    dplyr::arrange(dplyr::desc(dplyr::row_number()))
  
}

#' Process NGT records
#'
#' Process the NGT and associated oxygen isotope records, which includes reading
#' in the data, skipping unsuitabe records, and calculating anomaly time series.
#'
#' @param path file path of the NGT data set. For \code{NULL} (the default) the
#'   default file path set by the reading function \code{readNGT} is used.
#' @param reference.period numeric vector of ages of the reference period for
#'   calculating the anomaly time series; defaults to the standard 1961-1990
#'   period.
#' @return a data frame of anomaly time series of the relevant NGT and
#'   associated isotope records.
#' @author Thomas Münch
#'
processNGT <- function(path = NULL, reference.period = 1990 : 1961) {

  # Read data
  if (is.null(path)) {
    ngt <- readNGT()
  } else {
    ngt <- readNGT(path = path)
  }

  # Produce anomaly time series
  ngt <- makeAnomalies(ngt, reference.period = reference.period)

  return(ngt)

}

#' Process NGT accumulation data
#'
#' Process the NGT and associated accumulation rate records, which includes
#' reading in the data and calculating anomaly time series.
#'
#' @param path file path of the accumulation data set relative to the 'exNGT'
#'   repository.
#' @param reference.period numeric vector of ages of the reference period for
#'   calculating the anomaly time series; defaults to the standard 1961-1990
#'   period.
#' @return a data frame of the accumulation rate anomaly time series.
#' @author Thomas Münch
#'
processAccumulation <- function(path = "data/NGT2012_Accumulation_20220520.csv",
                                reference.period = 1990 : 1961) {

  if (file.exists(path)) {
    ngt <- read.csv(path, header = TRUE)
  } else {
    stop("No such file or directory; check path or permissions.")
  }

  ngt <- makeAnomalies(ngt, reference.period = reference.period)

  return(ngt)

}

#' Read DMI instrumental temperatures
#'
#' Read the Danish Meteorological Insitute (DMI) instrumental temperature
#' records of the Greenlandic stations Pituffik, Upernavik and Danmarkshavn.
#'
#' @param path file path (relative to working directory) of the DMI data
#'   set.
#' @return a data frame of four columns and 146 rows with the read DMI
#'   data set for the three stations.
#' @author Thomas Münch
#'
readDMI <- function(path = "data/gr_annual_temperature_1873_2015.csv") {

  colNames <- c("Year", "Pituffik", "Upernavik", "Danmarkshavn")

  dat <- read.csv(path, header = TRUE, na.strings = "#NULL!")
  dat <- dat[nrow(dat) : 1, c(1, 2, 4, 12)]

  colnames(dat) <- colNames

  return(dat)

}

#' Read HadCrut instrumental temperature for the Arctic region
#'
#' Read the annual mean HadCrut v5.0.1 60 to 90 degree N area mean temperature
#' data set from 1850 to 2021 CE.
#'
#' @param path file path (relative to working directory) of the HadCrut data
#'   set.
#' @return a data frame of two columns and 172 rows with the read HadCrut
#'   data set.
#' @author Thomas Münch
#'
readHadCrut <- function(path = "data/HadCrut_5.0.1_Arctic_annual.csv") {

  read.csv(path, header = TRUE) %>%
    dplyr::arrange(dplyr::desc(dplyr::row_number()))

}

#' Extend Arctic2k by HadCrut
#'
#' Extend the Arctic2k annual mean temperature reconstruction data set until
#' 2011 CE using the annual mean HadCrut v5.0.1 60 to 90 degree N area mean
#' temperature data.
#'
#' @param a2k data frame with the read Arctic2k reconstruction data obtained
#'   from \code{readArctic2k()}.
#' @param path file path (relative to working directory) of the HadCrut data
#'   set.
#' @return a data frame of four columns and 2011 rows with the Arctic2k
#'   data set extended by HadCrut until 2011 CE. Note that the columns
#'   \code{"2SigmaLow"} and \code{"2SigmaHigh"} are \code{NA} for the extension
#'   period.
#' @author Thomas Münch
#'
extendWithHadCrut <- function(a2k,
                              path = "data/HadCrut_5.0.1_Arctic_annual.csv") {

  # ensure identical anomalies 1961-1990
  d <- mean(subsetData(readHadCrut(), 1990 : 1961, "TempAnomaly")) -
    mean(subsetData(a2k, 1990 : 1961, "TempAnomaly"))

  readHadCrut() %>%
    dplyr::mutate(TempAnomaly = TempAnomaly - d) %>%
    dplyr::mutate(`2SigmaLow` = NA, `2SigmaHigh` = NA) %>%
    dplyr::filter(Year >= 2001 & Year <= 2011) %>%
    rbind(a2k)

}

#' Read MAR3.5 temperature and melt
#'
#' Read the MAR3.5 2m monthly surface air temperature averaged over the
#' NGT-2012 firn core locations as well as the Greenland melt rate, and process
#' the data to annual mean anomalies relative to the 1961-1990 CE reference
#' interval.
#'
#' @param path file path (relative to working directory) of the MAR data set.
#' @return a data frame of three columns and 141 rows with the read MAR3.5
#'   temperature and runoff (inverse melt rate) data set from 1871 to 2011 CE.
#' @author Thomas Münch
#'
readMAR <- function(path = "data/MAR_3.5_t2m_melt_annual.csv") {

  read.csv(path, header = TRUE) %>%
    setNames(c("Year", "t2m", "melt")) %>%
    dplyr::arrange(dplyr::desc(Year)) %>%
    dplyr::mutate(melt = -1 * melt) %>% # convert to runoff
    as.data.frame() %>%
    makeAnomalies()

}

#' Read GBI data
#'
#' Read the Greenland Blocking Index (GBI) data of Hanna et al. (2016) for the
#' 1851-2011 CE time period.
#'
#' @param path file path (relative to working directory) of the GBI data set.
#' @return a data frame of six columns and 161 rows with the read GBI annual and
#'   seasonal data from 1851 to 2011 CE.
#' @author Thomas Münch
#'
readGBI <- function(path = "data/gbi.monthly.csv") {

  read.csv(file = path, skip = 4) %>%
    dplyr::arrange(dplyr::desc(Year))

}

#' Read 20CR Greenland data
#'
#' Read the annual mean, area-weighted mean NOAA-CIRES-DOE 20CR v3 reanalysis
#' temperature data set from 1836 to 2015 CE for the Greenland region of the NGT
#' ice cores (72 to 80 degree latitude, -51 to -36 longitude).
#'
#' @param path file path (relative to working directory) of the 20CR data set.
#' @return a data frame of two columns and 180 rows with the read 20CRv3
#'   Greenland temperature time series from 1836 to 2015 CE.
#' @author Thomas Münch
readTwenCR <- function(path = "data/NOAA_CIRES_DOE_20CR_v3_annual_ngtregion_areamean.csv") {

  read.csv(file = path, header = TRUE) %>%
    dplyr::arrange(dplyr::desc(Year))

}

#' Site climatological parameters
#'
#' This function provides the relevant climatological parameters at the firn
#' core sites to estimate density and diffusion length profiles.
#'
#' @return A tibble with 16 rows and 6 columns. The columns list for each firn
#'   core site in this order the
#'   \itemize{
#'     \item number,
#'     \item site name,
#'     \item elevation in m a.s.l.,
#'     \item mean annual temperature in degree Celsius,
#'     \item mean annual accumulation rate in kg/m^2/yr,
#'     \item average surface snow density in kg/m^3, and
#'     \item atmospheric surface pressure in mbar.
#'   }
#'   The surface pressure is calculated from elevation and temperature using the
#'   barometric formula.
#'
#' @references
#'
#' * B16-B30 elevation, temperature and accumulation rate data from
#' Weissbach et al., Clim. Past, https://doi.org/10.5194/cp-12-171-2016, 2016.
#'
#' * GRIP and NGRIP elevation, temperature and accumulation rate data from
#' Vinther et al., J. Geophys. Res., https://doi.org/10.1029/2005JD006921, 2006.
#'
#' * GISP2 elevation and temperature data from
#' Grootes et al., Nature, https://doi.org/10.1038/366552a0, 1993.
#'
#' * GISP2 accumulion rate data estimated from Fig. 4 in
#' Cuffey and Clow, J. Geophys. Res., https://doi.org/10.1029/96JC03981, 1997.
#'
#' * NEEM elevation, temperature and accumulation rate data from
#' NEEM community members, Nature, https://doi.org/10.1038/nature11789, 2013.
#'
#' * NEGIS elevation and accumulation rate data from
#' Vallelonga et al., The Cryosphere, https://doi.org/10.5194/tc-8-1275-2014, 2014.
#'
#' * NEGIS temperature from
#' Zuhr et al., The Cryosphere Discuss., https://doi.org/10.5194/tc-2021-36, 2021.
#'
#' * Average surface snow density (0-1 m) from
#' Schaller et al., PANGAEA, https://doi.org/10.1594/PANGAEA.867874, 2016.
#' @author Thomas Münch
#'
loadClimPar <- function() {

  barometricFormula <- function (T, h, p0 = 1013) {

    kR <- 8.314
    kg <- 9.807
    kM <- 0.02896

    T  <- T + 273.15
    hs <- kR * T / (kg * kM)

    p0 * exp(-h / hs)
  }

  climPar <- tibble::tribble(
    ~Line, ~`Site`, ~`Elevation`, ~`meanTemperature`, ~`accRate`, ~`surfaceDensity`,
    # -- / ------ / ----------- / ----------------- / --------- / ---------------- /
    1,     "B16",   3040,         -32.5,              141,        310,
    2,     "B17",   2820,         -32.3,              114,        310,
    3,     "B18",   2508,         -32.3,              103,        310,
    4,     "B20",   2147,         -30.9,              98,         310,
    5,     "B21",   2185,         -30.1,              109,        310,
    6,     "B22",   2242,         -29.8,              145,        310,
    7,     "B23",   2543,         -29.3,              121,        310,
    8,     "B26",   2598,         -30.3,              176,        310,
    9,     "B27.28",2733,         -30.6,              180,        310,
    10,    "B29",   2874,         -31.6,              149,        310,
    11,    "B30",   2947,         -31.8,              166,        310,
    12,    "NGRIP", 2917,         -32,                175,        310,
    13,    "GISP2", 3208,         -31,                220,        310,
    14,    "GRIP",  3230,         -32,                212,        310,
    15,    "NEGIS", 2702,         -29,                101,        310,
    16,    "NEEM",  2450,         -29,                202,        310
  )

  dplyr::mutate(climPar,
                surfacePressure = barometricFormula(meanTemperature, Elevation))

}

#' Site coordinates
#'
#' This function provides the coordinates of the used firn core and weather
#' station sites.
#'
#' @return  A tibble with 19 rows and 4 columns. For each site, the columns list
#'   in this order the
#'   \itemize{
#'     \item number,
#'     \item site name,
#'     \item latitude in degree North, and
#'     \item longitude in degree East.
#' @author Thomas Münch
#'
loadPositions <- function() {

  tibble::tribble(
    ~Line, ~Identifier, ~`Site`, ~`Latitude`, ~`Longitude`,
    # -- / ---------- / ------ / ---------- / ----------- /
    1,     "core",      "B16",   73.94,       -37.63,
    2,     "core",      "B17",   75.25,       -37.63,
    3,     "core",      "B18",   76.62,       -36.40,
    4,     "core",      "B20",   78.83,       -36.50,
    5,     "core",      "B21",   80.00,       -41.14,
    6,     "core",      "B22",   79.34,       -45.91,
    7,     "core",      "B23",   78.00,       -44.00,
    8,     "core",      "B26",   77.25,       -49.22,
    9,     "core",      "B27/28",76.66,       -46.82,
    10,    "core",      "B29",   76.00,       -43.50,
    11,    "core",       "B30",   75.00,       -42.00,
    12,    "core",       "NGRIP", 75.10,       -42.30,
    13,    "core",       "GISP2", 72.60,       -38.50,
    14,    "core",       "GRIP",  72.60,       -37.60,
    15,    "core",       "NEGIS", 75.62,       -35.96,
    16,    "core",       "NEEM",  77.45,       -51.06,
    17,    "station",    "Pituffik",     76.53,   -68.75,
    18,    "station",    "Upernavik",    72.78,   -56.15,
    19,    "station",    "Danmarkshavn", 76.77,   -18.67
    )

}

#' Count records
#'
#' Count the number of available NGT records per year.
#'
#' @param ngt data frame with the NGT isotope data.
#' @return data frame with the NGT age scale as the first and the number of
#'   available isotope records as the second column.
#' @author Thomas Münch
#'
countRecords <- function(ngt) {

  data.frame(Year = ngt$Year,
             n = apply(ngt, MARGIN = 1, function(x) sum(!is.na(x))) - 1)

}

#' Anomaly time series
#'
#' Calculate anomaly time series relative to the mean of a given reference
#' period.
#'
#' @param data a data set for which to calculate the anomalies; can be a single
#'   numeric vector or a data frame or matrix of vectors. In the latter case,
#'   the first column can be used to provide the age scale of the data; else use
#'   the \code{age} parameter for this.
#' @param age provide the age scale of \code{data} here if it is not supplied
#'   in the latter parameter itself.
#' @param reference.period numeric vector with the ages of the reference period.
#' @param hasAgeColumn logical; do the supplied \code{data} contain an age
#'   column as the first column? Setting this to \code{FALSE} but leaving
#'   \code{age} as \code{NULL} will cause an error.
#' @return numeric vector, matrix or data frame with the anomaly time series.
#' @author Thomas Münch
#'
makeAnomalies <- function(data, age = NULL, reference.period = 1990 : 1961,
                          hasAgeColumn = TRUE) {

  if (is.null(age) & !hasAgeColumn) {
    stop("Age scale information missing.", call. = FALSE)
  }

  if (hasAgeColumn) {
    if (is.null(dim(data))) {
      stop("'data' is not a matrix or data frame.", call. = FALSE)
    }
    if (ncol(data) < 2) {
      stop("'data' only has 1 column; need at least 2.", call. = FALSE)
    }
    colNames <- colnames(data)
    Year <- data[, 1]
    data <- data[, -1]
  } else {
    Year <- age
  }

  i <- match(reference.period, Year)

  if (is.null(dim(data))) {

    result <- data - mean(data[i], na.rm = TRUE)

  } else {

    result <- apply(data, 2, function(x) {
      x - mean(x[i], na.rm = TRUE)})

  }

  if (hasAgeColumn) {
    result <- data.frame(Year, result)
    colnames(result) <- colNames

  }

  return(result)
}

#' Subset data set
#'
#' Subset a data set by specifying a time window and variable name.
#'
#' @param x a data frame or tibble including a time column and one or more named
#'   data columns.
#' @param t vector of time points to subset from \code{x}.
#' @param var the name of the data column to extract from \code{x}.
#' @param timeColumn the name of the time column in \code{x}; defaults to
#'   "Year".
#' @return a numeric vector with the data values for the requested data variable
#'   and the specified time window, ordered in time as in the original data
#'   set \code{x}.
#' @author Thomas Münch
#' @examples
#' NGT <- processNGT()
#' subsetData(NGT, t = 2011 : 2000, var = "B18_12")
#'
subsetData <- function(x, t, var, timeColumn = "Year") {

  x %>%
    dplyr::filter(!!as.name(timeColumn) %in% t) %>%
    dplyr::pull(!!as.name(var))
}

#' Isotope data for spectral analyses
#'
#' This is a wrapper function to quickly select those North Greenland isotope
#' records which have continuous data in a given time window; the default time
#' window is used to select those records suitable for the applied spectral
#' analyses.
#'
#' @param timeWindow vector of time points; default range is the one suited for
#'   the spectral analyses applied in the paper.
#' @return a data frame with all the North Greenland isotope data from the data
#'   compilation that have continuous data in the specified time window.
#' @author Thomas Münch
#'
selectNGTForSpectra <- function(timeWindow = 1979 : 1505) {

  noMissingVal <- function(x) {!any(is.na(x))}

  processNGT() %>%
    stackNGT(stack = FALSE) %>%
    dplyr::mutate(B23 = approx(Year, B23, Year)$y) %>%
    dplyr::slice(match(timeWindow, Year)) %>%
    dplyr::select(where(noMissingVal)) %>%
    dplyr::select(-Year)
}

#' Compile data for histogram
#'
#' This wrapper function compiles the NGT-2012 stacked data for the histogram
#' analysis depending on the data source and chosen stacking and merging
#' methods.
#'
#' @param data character string to signal the data source: one of "iso" (d18O
#'   isotope data) or "acc" (accumulation data).
#' @param type character string to signal the stacking method: one of "main"
#'   (main paper analyses), "stack_old_new" or "stack_all".
#' @param filter.window single integer giving the size of the running mean
#'   window to use for filtering the stack data.
#' @param adjustMean logical; whether to adjust the mean upon merging; see the
#'   description in \code[{mergeCores()}.
#' @param mergePoint character signalling the time point of merging; see the
#'   description in \code{mergeCores()}.
#' @param use_NEGIS_NEEM logical; whether to include the NEGIS and NEEM records
#'   in the stack.
#' @param diffuse logical; whether to forward diffuse the NGT data prior to
#'   the merging and stacking in order to mimic maximum smoothing at the firn
#'   ice transition.
#' @param nfix non-negative integer to set a constant number of records which
#'   are selected from the merged NGT-2012 data for stacking ("frozen" stack).
#' @return a data frame the running mean filtered stack data.
#' @author Thomas Münch
#'
selectHistogramData <- function(data = "iso", type = "main", filter.window = 11,
                                adjustMean = TRUE, mergePoint = "start",
                                use_NEGIS_NEEM = TRUE, diffuse = FALSE,
                                nfix = NULL) {

  if (!data %in% c("iso", "acc")) {
    stop("'data' must be one of 'iso' or 'acc'.", call. = FALSE)
  }

  if (!type %in% c("main", "stack_old_new", "stack_all", "fix_N")) {
    stop("'type' must be one of 'main', 'stack_old_new', 'stack_all'",
         " or 'fix_N'.", call. = FALSE)
  }

  if (length(nfix) & type != "fix_N") {
    stop("Fixed number of records incompatible with",
         " this stacking method, use 'fix_N'.", call. = FALSE)
  }

  if (!length(nfix) & type == "fix_N") {
    stop("Method 'fix_N' requires setting a fixed number of records.",
         call. = FALSE)
  }

  if (data == "acc" & type != "stack_all") {
    stop("Only type 'stack_all' is supported for accumulation data.",
         call. = FALSE)
  }

  selectFixedNumber <- function(NGT, nfix) {

    NGT <- NGT %>%
      dplyr::mutate(B18 = approx(Year, B18, Year)$y)

    rslt <- NULL
    nrow <- nrow(NGT)

    for (i in 1 : nrow) {

      row <- c(na.omit(unlist(unname(NGT[i, ]))))

      if (length(row) >= (nfix + 1)) {
        rslt <- rbind(rslt, row[1 : (nfix + 1)])
      } else {
        break
      }
    }

    return(rslt)
  }

  if (data == "iso") {
    NGT <- processNGT()
  } else {
    NGT <- processAccumulation()
  }

  if (diffuse) {

    cat("Forward diffusing data...\n")
    NGT <- diffuseNGT(NGT)
    cat("done.\n")
  }

  if (length(nfix)) {

    NGT <- NGT %>%
      filterData(window = filter.window) %>%
      stackNGT(stack = FALSE) %>%
      selectFixedNumber(nfix = nfix)
  }

  switch(type,

         main = filterData(NGT, window = filter.window) %>%
           mergeCores(adjustMean = adjustMean, mergePoint = mergePoint) %>%
           stackExtendedCores(filterData(NGT, window = filter.window)),

         stack_old_new = NGT %>%
           stackOldAndNew(use_NEGIS_NEEM = use_NEGIS_NEEM) %>%
           filterData(window = filter.window) %>%
           mergeCores(sites = "stack", adjustMean = adjustMean,
                      mergePoint = mergePoint),

         stack_all = NGT %>%
           filterData(window = filter.window) %>%
           stackAllCores(use_NEGIS_NEEM = use_NEGIS_NEEM),

         fix_N = data.frame(Year = NGT[, 1], stack = rowMeans(NGT[, -1]))
         )
}
