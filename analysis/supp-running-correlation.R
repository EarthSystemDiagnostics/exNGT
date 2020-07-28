##
## aim:
## analyse running correlation between the NGT stack and Arctic2k temperature
## relation:
## NGT paper supplementary; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

library(dplyr)

# ------------------------------------------------------------------------------
# Processing parameters

filter.window <- 11
permil2temperature <- 1 / 0.67
correlation.window <- 101

analysis.period <- c(1000, 2000)

# ------------------------------------------------------------------------------
# Load NGT and Arctic2k data and process for correlation analysis

NGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT() %>%
  dplyr::mutate(stack = permil2temperature * stack) %>%
  dplyr::filter(Year >= analysis.period[1] & Year <= analysis.period[2])

Arctic2k <- readArctic2k() %>%
  filterData(window = filter.window) %>%
  dplyr::filter(Year >= analysis.period[1] & Year <= analysis.period[2])

# ------------------------------------------------------------------------------
# Do correlation analysis

correlation <- doRunningCorrelation(NGT$stack, Arctic2k$TempAnomaly,
                                    window = correlation.window)

# ------------------------------------------------------------------------------
# Plot

ylim <- c(-1, 1)

xlab <- "Year CE"
ylab <- "Correlation"

time <- analysis.period[1] : analysis.period[2]

Quartz(file = "./fig/ngt-arctic2k-running-correlation.pdf", height = 4.5)

plot(time, correlation, type = "n",
     axes = FALSE, xlab = "", ylab = "", ylim = ylim)

axis(1)
axis(2)

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

lines(time, rep(0, length(time)), lty = 3, lwd = 1.5)

lines(time, correlation, col = "black", lwd = 2.5)

dev.off()
