##
## aim:
## analyse the histograms of all stacking possibilities.
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

# ------------------------------------------------------------------------------
# Calculate stacks

# filtered stack of old and new records
stack_old_new <- NGT %>%
  stackOldAndNew() %>%
  filterData(window = filter.window)

# filtered NGT data set
NGT.filtered <- NGT %>%
  filterData(window = filter.window)

stacks <- list(

  single_extension_naM_start = NGT.filtered %>%
    mergeCores(adjustMean = FALSE, method = 1) %>%
    stackExtendedCores(NGT.filtered),

  single_extension_naM_end = NGT.filtered %>%
    mergeCores(adjustMean = FALSE, method = 2) %>%
    stackExtendedCores(NGT.filtered),

  single_extension_aM_start = NGT.filtered %>%
    mergeCores(adjustMean = TRUE, method = 1) %>%
    stackExtendedCores(NGT.filtered),

  single_extension_aM_end = NGT.filtered %>%
    mergeCores(adjustMean = TRUE, method = 2) %>%
    stackExtendedCores(NGT.filtered),

  old_new_naM_start = stack_old_new %>%
    mergeCores(sites = "stack", adjustMean = FALSE, method = 1),

  old_new_naM_end = stack_old_new %>%
    mergeCores(sites = "stack", adjustMean = FALSE, method = 2),

  old_new_aM_start = stack_old_new %>%
    mergeCores(sites = "stack", adjustMean = TRUE, method = 1),

  old_new_aM_end = stack_old_new %>%
    mergeCores(sites = "stack", adjustMean = TRUE, method = 2),

  single_cores = NGT %>%
    stackAllCores() %>%
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
                               breaks = seq(-0.2, 0.2, 0.01),
                               xlim = c(-0.1, 0.1), ylim = c(0, 0.15),
                               xlab = "slope")

dev.off()

