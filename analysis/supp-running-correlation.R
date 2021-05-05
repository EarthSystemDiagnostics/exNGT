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
correlation.window <- 101

analysis.period <- c(1000, 2000)

# ------------------------------------------------------------------------------
# Load NGT and Arctic2k data and process for correlation analysis

NGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()

Arctic2k <- readArctic2k() %>%
  filterData(window = filter.window)

NGTsubset <- NGT %>%
  dplyr::filter(Year >= analysis.period[1] & Year <= analysis.period[2])

Arctic2kSubset <- Arctic2k %>%
  dplyr::filter(Year >= analysis.period[1] & Year <= analysis.period[2])

# ------------------------------------------------------------------------------
# Do correlation analysis

correlation <- doRunningCorrelation(NGTsubset$stack, Arctic2kSubset$TempAnomaly,
                                    window = correlation.window)

# ------------------------------------------------------------------------------
# Plot

time <- analysis.period[1] : analysis.period[2]

xlim <- range(time)
ylim <- c(-1, 1)

ylim.ngt <- c(-4, 2)
ylim.a2k <- c(-2, 4)

xlab <- "Year CE"
ylab <- "Correlation"

ylab.ngt <- bquote(delta^{"18"} * "O anomaly (\u2030)")
ylab.a2k <- bquote("Temperature anomaly" * " (" * degree * "C)")

x <- 2170
y <- -0.5

col <- c("black", "dodgerblue4")

grfxtools::Quartz(file = "./fig/ngt-arctic2k-running-correlation.pdf",
                  height = 7, width = 8.9, mfrow = c(2, 1),
                  mar = c(0, 0, 0, 0), oma = c(5, 5, 0.5, 5))

plot(NGT, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim.ngt)

axis(2, at = seq(-2, 2, 1))
mtext(ylab.ngt, side = 2, line = 3.25,
      cex = par()$cex.lab * par()$cex, las = 0, adj = 1)

mtext("a", side = 3, adj = 0.01, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab, col = col[1])
mtext("b", side = 3, adj = 0.99, padj = 8,
      line = -1, font = 2, cex = par()$cex.lab, col = col[2])

abline(h = 0, lty = 3, lwd = 1.5, col = col[1])
lines(NGT, type = "l", col = col[1], lwd = 2.5)

par(new = TRUE)

plot(Arctic2k$Year, Arctic2k$TempAnomaly, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = xlim, ylim = ylim.a2k)

axis(4, at = seq(-2, 1, 1), col = col[2], col.axis = col[2])
text(x, y, ylab.a2k, srt = -90, xpd = NA, cex = par()$cex.lab * par()$cex,
     col = col[2])

abline(h = 0, lty = 3, lwd = 1.5, col = col[2])
lines(Arctic2k$Year, Arctic2k$TempAnomaly, col = col[2], lwd = 2.5)

plot(time, correlation, type = "n",
     axes = FALSE, xlab = "", ylab = "", xlim = xlim, ylim = ylim)

axis(1)
axis(2)

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

mtext("c", side = 3, adj = 0.01, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab, col = col[1])

lines(time, rep(0, length(time)), lty = 3, lwd = 1.5)

lines(time, correlation, col = "black", lwd = 2.5)

dev.off()
