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

# ------------------------------------------------------------------------------
# Processing parameters

filter.window <- 11
correlation.window <- 101

analysis.period <- 2000 : 1000

# ------------------------------------------------------------------------------
# Load NGT and Arctic2k data and process for correlation analysis

NGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()

Arctic2k <- readArctic2k() %>%
  filterData(window = filter.window)

NGTsubset <- processNGT() %>%
  stackNGT() %>%
  subsetData(analysis.period, "stack")

Arctic2kSubset <- readArctic2k() %>%
  subsetData(analysis.period, "TempAnomaly")

# Number of records contributing to NGT-2012 stack
nRecords <- processNGT() %>%
  stackNGT(stack = FALSE) %>%
  countRecords()

# ------------------------------------------------------------------------------
# Do correlation analysis

correlation <- estimateRunningCorrelation(
  NGTsubset, Arctic2kSubset, nmc = 10000,
  correlation.window = correlation.window, filter.window = filter.window)

# ------------------------------------------------------------------------------
# Plot

xlim <- range(analysis.period)
ylim <- c(-0.5, 1)

ylim.ngt <- c(-4, 2)
ylim.nbr <- c(-130, 20)
ylim.a2k <- c(-2, 4)

xlab <- "Year CE"
ylab <- "Correlation"

ylab.ngt <- grfxtools::LabelAxis("NGT-2012")
ylab.a2k <- grfxtools::LabelAxis("Arctic2k", unit = "celsius")

x1 <- 845
x2 <- 2160
y1 <- 0.
y2 <- 10
y3 <- -0.5

col <- c("black", "dodgerblue4", "darkgrey")

grfxtools::Quartz(file = "./fig/supplement-ngt-arctic2k-running-correlation.pdf",
                  height = 7, width = 8.9, mfrow = c(2, 1),
                  mar = c(0, 0, 0, 0), oma = c(5, 5, 0.5, 5))

plot(NGT, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim.ngt)

axis(2, at = seq(-2, 2, 1))
text(x1, y1, ylab.ngt, srt = +90, xpd = NA, cex = par()$cex.lab * par()$cex,
     col = col[1])

mtext("a", side = 3, adj = 0.01, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab, col = col[1])
mtext("b", side = 3, adj = 0.99, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab, col = col[3])
mtext("c", side = 3, adj = 0.99, padj = 8,
      line = -1, font = 2, cex = par()$cex.lab, col = col[2])

abline(h = 0, lty = 2, lwd = 1.5, col = "darkgrey")
lines(NGT, type = "l", col = col[1], lwd = 2.5)

par(new = TRUE)

plot(nRecords, type = "l", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim.nbr, col = col[3], lwd = 2)

axis(4, at = c(0, 10, 20), col = col[3], col.axis = col[3])
text(x2, 10, "N", srt = -90, xpd = NA, cex = par()$cex.lab,
     col = col[3])

par(new = TRUE)

plot(Arctic2k$Year, Arctic2k$TempAnomaly, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = xlim, ylim = ylim.a2k)

axis(4, at = seq(-2, 1, 1), col = col[2], col.axis = col[2])
text(x2, y3, ylab.a2k, srt = -90, xpd = NA, cex = par()$cex.lab * par()$cex,
     col = col[2])

abline(h = 0, lty = 2, lwd = 1.5, col = "darkgrey")
lines(Arctic2k$Year, Arctic2k$TempAnomaly, col = col[2], lwd = 2.5)

plot(analysis.period, correlation$dat, type = "n",
     axes = FALSE, xlab = "", ylab = "", xlim = xlim, ylim = ylim)

axis(1)
axis(2)

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

mtext("d", side = 3, adj = 0.01, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab, col = col[1])

lines(analysis.period, rep(0, length(analysis.period)),
      lty = 2, lwd = 1.5, col = "darkgrey")

lines(analysis.period, correlation$dat, col = "black", lwd = 2.5)
lines(analysis.period, correlation$r.upper, lty = 2, lwd = 1, col = "black")
lines(analysis.period, correlation$r.lower, lty = 2, lwd = 1, col = "black")

dev.off()

save(correlation,
     file = "./out/supplement-ngt-arctic2k-running-correlation.rda")
