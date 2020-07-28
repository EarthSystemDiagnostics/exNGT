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
plot_sub_periods <- FALSE

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

xanml <- c(800, 2020)
yanml <- rep(0, 2)

startNew <- 1993
i <- match(startNew, stackedNGT$Year)
n <- nrow(stackedNGT)

periodColors = c("#e41a1c", "#ff7f00", "#377eb8", "#4daf4a", "#984ea3")

periods <- list(
  p1 = 2011 : 1997, p2 = 1996 : 1982,
  p3 = 1938 : 1924, p4 = 1884 : 1870,
  p5 = 1424 : 1410
)


if (plot_record_number) {

  # for plot with number of records
  Quartz(file = "./fig/ngt-record-number-time-periods.pdf",
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

if (plot_sub_periods) {

  for (iT in 1 : length(periods)) {

    nt <- length(periods[[iT]])
    polygon(x = c(periods[[iT]], rev(periods[[iT]])),
            y = c(rep(-3, nt), rep(2.5, nt)),
            col = adjustcolor(periodColors[iT], alpha = 0.4),
            border = NA)

  }
}

lines(xanml, yanml, lty = 3, lwd = 1.5)

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

lines(xanml, yanml, lty = 3, lwd = 1.5)

lines(stackedTemperatureNGT, col = "darkgrey")

lines(filteredStackedTemperatureNGT, col = col[1], lwd = 2.5)

lines(Arctic2k$Year, Arctic2k$TempAnomaly,
      col = adjustcolor(col[2], alpha = 0.6))
lines(filteredArctic2k$Year, filteredArctic2k$TempAnomaly, col = col[2], lwd = 2.5)

legend("topleft", c("NGT stack", "Arctic2k"), col = col, lwd = 2.5, bty = "n")

dev.off()

# ------------------------------------------------------------------------------
# Fig. 3, anomaly and slope histogram for main stack

lab0 <- "Full data"
lab1 <- "1997 to 2011"
lab2 <- "1982 to 1996"
lab3 <- "1924 to 1938"
lab4 <- "1870 to 1884"
lab5 <- "1410 to 1424"

col = c("#e41a1c", "#ff7f00", "#377eb8", "#4daf4a", "#984ea3")

slopesFilteredStackedNGT <- estimateSlopes(filteredStackedNGT)

Quartz(file = "./fig/main-05-histograms.pdf",
       height = 5, width = 12)

layout(matrix(1 : 3, 1, 3), widths = c(0.43, 0.43, 0.14))
par(cex = 1)

plotHistogram(filteredStackedNGT, plot.legend = FALSE)
mtext("a", side = 3, line = -0.5, cex = par()$cex.lab * par()$cex,
      las = 1, font = 2, adj = 0.02, padj = 0.5)

plotHistogram(slopesFilteredStackedNGT,
              breaks = seq(-0.2, 0.2, 0.01), range = "pos",
              xlim = c(-0.1, 0.1), ylim = c(0, 0.2),
              xlab = "slope", plot.legend = FALSE)
mtext("b", side = 3, line = -0.5, cex = par()$cex.lab * par()$cex,
      las = 1, font = 2, adj = 0.02, padj = 0.5)

par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("topleft", c(lab0, lab1, lab2, lab3, lab4, lab5),
       col = c(adjustcolor(1, 0.2), adjustcolor(col, 0.6)),
       lty = 1, lwd = 10, bty = "n")

dev.off()
