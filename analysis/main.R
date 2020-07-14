##
## aim:
## make analyses and plots for the main paper part.
## relation:
## NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

library(magrittr)

# ------------------------------------------------------------------------------
# Load NGT and Arctic2k data

NGT <- processNGT()
Arctic2k <- readArctic2k()

# ------------------------------------------------------------------------------
# Processing parameters

filter.window <- 11
permil2temperature <- 1 / 0.67

# ------------------------------------------------------------------------------
# Fig. 1, main, full NGT

plot_record_number <- FALSE

# stack NGT
stackedNGT <- NGT %>%
  stackNGT()

# filter and stack NGT
filteredStackedNGT <- NGT %>%
  filterData(window = filter.window) %>%
  stackNGT()

# number of records contributing to the stack
nRecords <- NGT %>%
  stackNGT(stack = FALSE) %>%
  countRecords()

xlab  <- "Year CE"
ylab <- bquote(delta^{"18"} * "O anomaly (\u2030)")

xlim <- c(1000, 2020)
ylim <- c(-2.5, 2.5)

startNew <- 1993
i <- match(startNew, stackedNGT$Year)
n <- nrow(stackedNGT)

if (plot_record_number) {

  # for plot with number of records
  Quartz(file = "./fig/main-01-ngt-with-record-number.pdf",
         height = 4.5 , width = 8.9)
  par(mar = c(5, 5, 0.5, 5))

} else {

  Quartz(file = "./fig/main-01-ngt.pdf", height = 4.5)
}

plot(stackedNGT, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim)

axis(1)
axis(2)

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

lines(stackedNGT, col = "darkgrey")

lines(filteredStackedNGT[i : n, ], col = "black", lwd = 2.5)
lines(filteredStackedNGT[1 : i, ], col = "firebrick3", lwd = 2.5)

if (plot_record_number) {

  par(new = TRUE)

  plot(nRecords, type = "l", axes = FALSE, xlab = "", ylab = "",
       xlim = xlim, ylim = c(0, 150), col = "darkgrey", lwd = 1.5)

  axis(4, at = c(0, 10, 20))
  text(2180, 10, "Number", srt = -90, xpd = NA, cex = par()$cex.lab)

}

dev.off()

# ------------------------------------------------------------------------------
# Fig. 1, minor, NGT and Arctic2k

# convert NGT stack to temperature
stackedTemperatureNGT <- stackedNGT %>%
  dplyr::mutate(stack = permil2temperature * stack)
filteredStackedTemperatureNGT <- filteredStackedNGT %>%
  dplyr::mutate(stack = permil2temperature * stack)

# filter Arctic2k
filteredArctic2k <- Arctic2k %>%
  filterData(window = filter.window)

col <- c("black", "dodgerblue4")

xlab <- "Year CE"
ylab <- bquote("Temperature anomaly" * " (" * degree * "C)")

xlim <- c(1840, 2020)
ylim <- c(-4, 4)

Quartz(file = "./fig/main-02-ngt-arctic2k-comparison.pdf", height = 4.5)

plot(stackedTemperatureNGT, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim)

axis(1, at = seq(xlim[1], xlim[2], 20))
axis(2, at = seq(ylim[1], ylim[2], 1))

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

lines(stackedTemperatureNGT, col = "darkgrey")

lines(filteredStackedTemperatureNGT, col = col[1], lwd = 2.5)

lines(Arctic2k$Year, Arctic2k$TempAnomaly,
      col = adjustcolor(col[2], alpha = 0.6))
lines(filteredArctic2k$Year, filteredArctic2k$TempAnomaly, col = col[2], lwd = 2.5)

legend("topleft", c("NGT stack", "Arctic2k"), col = col, lwd = 2.5, bty = "n")

dev.off()

