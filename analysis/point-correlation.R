#!/usr/bin/env Rscript

##
## aim:
## analyse point correlation between 20CR reanalysis field and NGT-2012/Arctic2k
## relation:
## NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##

# make script to also run on the terminal
if (interactive()) {

  # interactive usage
  path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
  path <- "~/programming/R/exNGT" #Thomas

  setwd(path)
  source("init.R")

} else {

  # script usage
  pwd <- Sys.getenv()[["PWD"]]
  setwd(pwd)
  setwd("../")
  cat("Running program in", getwd(), "\n\n")
  source("init.R")
  cat("\n\n")

}

# number of Monte Carlo iterations
nmc <- 1000

# ------------------------------------------------------------------------------
# Load data and filter

filter.window <- 11

NGT <- processNGT() %>%
  stackNGT()

filteredNGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()

Arctic2k <- readArctic2k() %>%
  dplyr::select("Year", "TempAnomaly")

filteredArctic2k <- readArctic2k() %>%
  filterData(window = filter.window) %>%
  dplyr::select("Year", "TempAnomaly")

load(file = "data/NOAA_CIRES_DOE_20CR_v3_50to90_annual_field.rda",
     verbose = TRUE)

filteredTwenCR <- TwenCR
filteredTwenCR$dat <- apply(filteredTwenCR$dat, 2, filterData,
                            window = filter.window,
                            hasAgeColumn = FALSE)

GreenlandTwenCR <- readTwenCR()
ArcticTwenCR <- readArcticTwenCR()

# ------------------------------------------------------------------------------
# Wrapper function

getPointCorrelation <- function(data, filteredData, fieldData, filter.window,
                                nmc, analysis.period) {

  run <- function(x, time, data, filteredData,
                  filter.window, nmc, analysis.period) {

    signal <- data.frame(Year = time, dat = x)

    estimateCorrelation(data, filteredData, signal,
                        filter.window = filter.window, nmc = nmc,
                        analysis.period = analysis.period) %>%
      simplify2array()

  }

  field <- fieldData$dat
  time  <- fieldData$time

  pointCor <- field %>%
    apply(2, run, time, data, filteredData,
          filter.window, nmc, analysis.period) %>%
    t() %>%
    as.data.frame()    

  pointCor$lat <- fieldData$lat
  pointCor$lon <- fieldData$lon

  pointCor %>%
    dplyr::relocate(lat, lon)

}

# ------------------------------------------------------------------------------
# Run analyses for 20CR field with NGT-2012 and Arctic2k time series

cat("\n")

cat(as.character(Sys.time()), "Running annual NGT-20CR...\n")
cor.ngt.ann <- getPointCorrelation(NGT, NGT, TwenCR,
                                   filter.window = 1, nmc = nmc,
                                   analysis.period = 1836 : 2000)
cat(as.character(Sys.time()), "Running annual A2k-20CR...\n")
cor.a2k.ann <- getPointCorrelation(Arctic2k, Arctic2k, TwenCR,
                                   filter.window = 1, nmc = nmc,
                                   analysis.period = 1836 : 2000)
cat(as.character(Sys.time()), "Running 11-yr NGT-20CR...\n")
cor.ngt.11y <- getPointCorrelation(NGT, filteredNGT, filteredTwenCR,
                                   filter.window = filter.window, nmc = nmc,
                                   analysis.period = 1836 : 2000)
cat(as.character(Sys.time()), "Running 11-yr A2k-20CR...\n")
cor.a2k.11y <- getPointCorrelation(Arctic2k, filteredArctic2k, filteredTwenCR,
                                   filter.window = filter.window, nmc = nmc,
                                   analysis.period = 1836 : 2000)

cat("\n")
cat(as.character(Sys.time()), "...done.", "\n\n")

cat("Saving data...\n\n")

list(ngt.ann = cor.ngt.ann, a2k.ann = cor.a2k.ann,
     ngt.11y = cor.ngt.11y, a2k.11y = cor.a2k.11y) %>%
  saveRDS(file = "out/point-correlations-20CR.rds")

# ------------------------------------------------------------------------------
# Run analyses for 20CR field with 20CR Greenland and Arctic aremean time series

cat("\n")

cat(as.character(Sys.time()), "Running Greenland correlation with field...\n")
cor.greenland <- getPointCorrelation(GreenlandTwenCR, GreenlandTwenCR, TwenCR,
                                     filter.window = 1, nmc = nmc,
                                     analysis.period = 1836 : 2015)
cat(as.character(Sys.time()), "Running Arctic correlation with field...\n")
cor.arctic <- getPointCorrelation(ArcticTwenCR, ArcticTwenCR, TwenCR,
                                  filter.window = 1, nmc = nmc,
                                  analysis.period = 1836 : 2015)

cat("\n")
cat(as.character(Sys.time()), "...done.", "\n\n")

cat("Saving data...\n\n")

list(greenland = cor.greenland, arctic = cor.arctic) %>%
  saveRDS(file = "out/point-correlations-20CR-with-20CR.rds")

cat("Good bye.\n\n")

