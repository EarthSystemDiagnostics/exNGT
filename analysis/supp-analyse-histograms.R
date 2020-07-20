##
## aim:
## analyse the histograms of all stacking possibilities and their sensitivity
## on different methods as well as on diffusional smoothing.
## relation:
## NGT paper supplementary; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

library(magrittr)

# ------------------------------------------------------------------------------
# Input

# all NGT anomaly records
NGT <- processNGT()

# set filter window
filter.window <- 11


# ==============================================================================
# Test different stacking options
# ==============================================================================

# ------------------------------------------------------------------------------
# Calculate stacks

# filtered stack of old and new records
stack_old_new <- NGT %>%
  stackOldAndNew(use_NEGIS_NEEM = FALSE) %>%
  filterData(window = filter.window)

# filtered NGT data set
NGT.filtered <- NGT %>%
  filterData(window = filter.window)

stacks <- list(

  single_extension_naM_start = NGT.filtered %>%
    mergeCores(adjustMean = FALSE, method = 1) %>%
    stackExtendedCores(NGT.filtered, use_NEGIS_NEEM = FALSE),

  single_extension_naM_end = NGT.filtered %>%
    mergeCores(adjustMean = FALSE, method = 2) %>%
    stackExtendedCores(NGT.filtered, use_NEGIS_NEEM = FALSE),

  single_extension_aM_start = NGT.filtered %>%
    mergeCores(adjustMean = TRUE, method = 1) %>%
    stackExtendedCores(NGT.filtered, use_NEGIS_NEEM = FALSE),

  single_extension_aM_end = NGT.filtered %>%
    mergeCores(adjustMean = TRUE, method = 2) %>%
    stackExtendedCores(NGT.filtered, use_NEGIS_NEEM = FALSE),

  old_new_naM_start = stack_old_new %>%
    mergeCores(sites = "stack", adjustMean = FALSE, method = 1),

  old_new_naM_end = stack_old_new %>%
    mergeCores(sites = "stack", adjustMean = FALSE, method = 2),

  old_new_aM_start = stack_old_new %>%
    mergeCores(sites = "stack", adjustMean = TRUE, method = 1),

  old_new_aM_end = stack_old_new %>%
    mergeCores(sites = "stack", adjustMean = TRUE, method = 2),

  single_cores = NGT %>%
    stackAllCores(use_NEGIS_NEEM = FALSE) %>%
    filterData(window = filter.window)
)

# ------------------------------------------------------------------------------
# Estimate the trend slopes for all stacks

slopes <- lapply(stacks, estimateSlopes)

# ------------------------------------------------------------------------------
# Plot the histograms

n <- length(stacks)

xmain = c("a. Non-adjusted mean, extend at start",
          "b. Non-adjusted mean, extend at end",
          "c. Adjusted mean, extend at start",
          "d. Adjusted mean, extend at end",
          "e. Non-adjusted mean, extend at start",
          "f. Non-adjusted mean, extend at end",
          "g. Adjusted mean, extend at start",
          "h. Adjusted mean, extend at end",
          "i.")
ymain = c("Extend single cores", rep(NA, 3),
          "Extend stacks", rep(NA, 3),
          "Average all cores")

Quartz(file = "./fig/histograms-anomalies.pdf", height = 12, width = 16)
par(mfrow = c(3, 4), mar = c(5, 9, 4, 0.5))

for (i in 1 : n) plotHistogram(stacks[[i]], xmain = xmain[i], ymain = ymain[i])

dev.off()

Quartz(file = "./fig/histograms-slopes.pdf", height = 12, width = 16)
par(mfrow = c(3, 4), mar = c(5, 9, 4, 0.5))

for (i in 1 : n) plotHistogram(slopes[[i]], xmain = xmain[i], ymain = ymain[i],
                               breaks = seq(-0.2, 0.2, 0.01), range = "pos",
                               xlim = c(-0.1, 0.1), ylim = c(0, 0.15),
                               xlab = "slope")

dev.off()


# ==============================================================================
# Test other method variants
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Analyse dependence of main stack histogram on filter method

# main filtered NGT stack
filteredStackedNGT <- NGT %>%
  filterData(window = filter.window) %>%
  stackNGT()

# do main stack for different filter end point constraints
filteredStackedNGT.method1 <- NGT %>%
  filterData(window = filter.window, method = 1) %>%
  stackNGT()
filteredStackedNGT.method3 <- NGT %>%
  filterData(window = filter.window, method = 3) %>%
  stackNGT()

# ------------------------------------------------------------------------------
# 2. Main stack for constant number of records (n = 6)

# get records that contribute to main stack
mainStackRecords <- NGT %>%
  filterData(window = filter.window) %>%
  stackNGT(stack = FALSE)

# select at least the number of records available for recent year
nmin <- 6
mainStackFixN <- mainStackRecords %>%
  apply(1, function(row) {row[which(!is.na(row))[1 : (nmin + 1)]]}) %>%
  t()
# <- these are always six records but their composition varies through time!

# stack this data set
mainStackFixN <- data.frame(Year = mainStackFixN[, 1],
                            stack = rowMeans(mainStackFixN[, -1]))

# ------------------------------------------------------------------------------
# 3. Main stack after full forward diffusion

# get diffusion length estimates in units of years for each site

library(FirnR)
climatePar <- loadClimPar()

sigma <- FirnR::TemporalDiffusionLength(nt = nrow(NGT),
                                        T = climatePar$meanTemperature + 273.15,
                                        P = climatePar$surfacePressure,
                                        bdot = climatePar$accRate,
                                        rho.surface = climatePar$surfaceDensity,
                                        names = climatePar$Site)

# get the maximum diffusion length in the ice
sigma.ice <- lapply(sigma, max) %>%
  as.data.frame()

# get the differential diffusion length between local (firn) and ice value
sigma.differential <- sigma
for (i in 2 : ncol(sigma)) {
  sigma.differential[, i] <- sqrt((sigma.ice[, i])^2 - (sigma[, i])^2)
}

applyDiffusion <- function(x, siteID, sigma) {

  siteID    <- strsplit(siteID, "_")[[1]][1]
  siteIndex <- match(siteID, colnames(sigma))
  inm       <- which(!is.na(x))

  print(siteID)
  print(siteIndex)

  x[inm] <- FirnR::DiffuseRecord(x[inm], sigma = sigma[seq(inm), siteIndex])

  return(x)

}

diffusedNGT <- NGT
for (i in 2 : ncol(NGT)) {

  diffusedNGT[, i] <- applyDiffusion(NGT[, i], siteID = colnames(NGT)[i],
                                     sigma = sigma.differential)
}

# build main stack from fully forward diffused records
diffusedFilteredStackedNGT <- diffusedNGT %>%
  filterData(window = filter.window) %>%
  stackNGT()

# ------------------------------------------------------------------------------
# Plot the histograms

xmain = c("a. NGT main stack",
          "b. NGT main stack; minimum norm",
          "c. NGT main stack; minimum roughness",
          "d. NGT main stack; fixed record number",
          "e. NGT main stack; full forward diffusion")

Quartz(file = "./fig/histograms-different-methods.pdf", height = 8, width = 12)
par(mfrow = c(2, 3), mar = c(5, 5, 4, 0.5))

plotHistogram(filteredStackedNGT, xmain = xmain[1])
plotHistogram(filteredStackedNGT.method1, xmain = xmain[2])
plotHistogram(filteredStackedNGT.method3, xmain = xmain[3])

plot.new()
plotHistogram(mainStackFixN, xmain = xmain[4])
plotHistogram(diffusedFilteredStackedNGT, xmain = xmain[5])

dev.off()

