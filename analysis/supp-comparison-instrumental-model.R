##
## aim:
## compare the NGT stack with Greenlandic instrumental and model data.
## relation:
## NGT paper supplementary; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

# ------------------------------------------------------------------------------
# Load data

filter.window <- 11

permil2temperature.mid <- 1 / 0.67 # Greenland spatial slope
permil2temperature.low <- 1 / 1.1  # Masson-Delmotte et al. (2015)
permil2temperature.hig <- 2.1      # Vinther et al. (2009)

# NGT-2012 stack
filteredStackedNGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()

filteredStackedTemperatureNGT.mid <- filteredStackedNGT %>%
  dplyr::mutate(stack = permil2temperature.mid * stack)
filteredStackedTemperatureNGT.low <- filteredStackedNGT %>%
  dplyr::mutate(stack = permil2temperature.low * stack)
filteredStackedTemperatureNGT.hig <- filteredStackedNGT %>%
  dplyr::mutate(stack = permil2temperature.hig * stack)

# DMI records
filteredDMI <- readDMI() %>%
  makeAnomalies() %>%
  filterData(window = filter.window) %>%
  dplyr::filter(Year >= 1901)

# MAR data
filteredMAR <- readMAR() %>%
  filterData(window = filter.window)

# Annual mean GBI data
filteredAnnualGBI <- readGBI() %>%
  filterData(window = filter.window) %>%
  dplyr::select("Year", gbi = "annual")

# ------------------------------------------------------------------------------
# Plot

grfxtools::Quartz(file = "./fig/supplement-instrumental-model-comparison.pdf",
                  height = 10.5 , width = 9, mar = c(5, 5, 0.5, 5))

## layout(matrix(c(1, 1, 1, 2 : 7), 3, 3, byrow = TRUE),
##        widths = c(1, 1, 1), heights = c(1.2, 0.9, 0.9))
layout(matrix(c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), 5, 5, 6, 6, 7, 7),
              3, 6, byrow = TRUE), widths = rep(1, 6), heights = rep(1, 3))
par(cex = 0.8)

# ------------------------------------------------------------------------------
# Part 1

filteredStackedNGT <- filteredStackedNGT %>%
  dplyr::filter(Year >= 1901)

col <- c("black", "#1b9e77", "#d95f02", "#7570b3")

xlab  <- "Year CE"
ylab1 <- grfxtools::LabelAxis("NGT-2012", suffix = "anomaly")
ylab2 <- grfxtools::LabelAxis("DMI anomaly", unit = "celsius")

xlim <- c(1900, 2020)
ylim1 <- c(-1, 1.5)
ylim2 <- c(-2, 3)

x <- 2050 
y <- mean(ylim2)

plot(filteredStackedNGT, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim1)

axis(1)
axis(2)

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

mtext("a", side = 3, adj = 0.02, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab * par()$cex)

lines(filteredStackedNGT, col = col[1], lwd = 2.5)

par(new = TRUE)

plot(filteredStackedNGT, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim2)

axis(4)

text(x, y, ylab2, srt = -90, xpd = NA, cex = par()$cex.lab)

for (i in 2 : 4) lines(filteredDMI$Year, filteredDMI[, i], col = col[i], lwd = 2.5)

legend("bottom", c("NGT-2012", colnames(filteredDMI)[-1]),
       cex = 0.95, col = col, lwd = 2.5, bty = "n", inset = c(0, -0.01))

# ------------------------------------------------------------------------------
# Part 2

ylab0 <- grfxtools::LabelAxis("Temperature anomaly", unit = "celsius")
ylab1 <- grfxtools::LabelAxis("NGT-2012 anomaly", unit = "celsius")
ylab2 <- grfxtools::LabelAxis("Runoff anomaly", unit = "Gt",
                              unit.type = "trend", time.unit = "yr")
ylab3 <- "GBI"

xlim1 <- c(1870, 2012)
xlim2 <- c(1850, 2012)
ylim1 <- c(-2.5, 3.5)
ylim2 <- c(-75, 330)
ylim3 <- c(-0.5, 3)

x1 <- 2052.5
x2 <- 2055
y1 <- 125
y2 <- 0.25

col <- c("black", "#7570b3", "#d95f02", "#1b9e77")

# plot 1 -----
plot(filteredStackedTemperatureNGT.mid, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = xlim1, ylim = ylim1)

axis(1)
axis(2)
mtext(xlab, side = 1, line = 3.25, cex = par()$cex.lab * par()$cex)
mtext(ylab0, side = 2, line = 2.75, cex = par()$cex.lab * par()$cex, las = 0)
mtext("b", side = 3, adj = 0.02, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab * par()$cex)

lines(filteredStackedTemperatureNGT.mid, lwd = 2, col = col[1])

grfxtools::Polyplot(filteredStackedTemperatureNGT.mid$Year,
                    filteredStackedTemperatureNGT.low$stack,
                    filteredStackedTemperatureNGT.hig$stack,
                    col = col[1], alpha = 0.2)
lines(filteredStackedTemperatureNGT.low, lwd = 1, col = col[1])
lines(filteredStackedTemperatureNGT.hig, lwd = 1, col = col[1])

lines(filteredMAR$Year, filteredMAR$t2m, col = col[2], lwd = 2)

legend("bottom", c("NGT-2012", "MAR3.5.2"), lty = 1,
       lwd = 2, col = col[1 : 2], bty = "n", inset = c(0, -0.01))

# plot 2 -----
plot(filteredStackedTemperatureNGT.mid, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = xlim1, ylim = ylim1)

axis(1)
axis(2)
mtext(xlab, side = 1, line = 3.25, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)
mtext("c", side = 3, adj = 0.02, padj = 0.5,
      line = -2, font = 2, cex = par()$cex.lab * par()$cex)

lines(filteredStackedTemperatureNGT.mid, lwd = 2, col = col[1])

grfxtools::Polyplot(filteredStackedTemperatureNGT.mid$Year,
                    filteredStackedTemperatureNGT.low$stack,
                    filteredStackedTemperatureNGT.hig$stack,
                    col = col[1], alpha = 0.2)
lines(filteredStackedTemperatureNGT.low, lwd = 1, col = col[1])
lines(filteredStackedTemperatureNGT.hig, lwd = 1, col = col[1])

par(new = TRUE)

plot(filteredMAR$Year, filteredMAR$melt, type = "l", axes = FALSE,
     col = col[3], lwd = 2, xlim = xlim1, ylim = ylim2,
     xlab = "", ylab = "")

axis(4, at = seq(-50, 300, 50), col = col[3], col.axis = col[3])
text(x1, y1, ylab2, srt = -90, xpd = NA, col = col[3], cex = par()$cex.lab)

# plot 3 -----
plot(filteredStackedTemperatureNGT.mid, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = xlim2, ylim = ylim1)

axis(1)
axis(2)
mtext(xlab, side = 1, line = 3.25, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 2.75, cex = par()$cex.lab * par()$cex, las = 0)
mtext("d", side = 3, adj = 0.02, padj = 0.5,
      line = -2, font = 2, cex = par()$cex.lab * par()$cex)

lines(filteredStackedTemperatureNGT.mid, lwd = 2, col = col[1])

grfxtools::Polyplot(filteredStackedTemperatureNGT.mid$Year,
                    filteredStackedTemperatureNGT.low$stack,
                    filteredStackedTemperatureNGT.hig$stack,
                    col = col[1], alpha = 0.2)
lines(filteredStackedTemperatureNGT.low, lwd = 1, col = col[1])
lines(filteredStackedTemperatureNGT.hig, lwd = 1, col = col[1])

par(new = TRUE)

plot(filteredAnnualGBI, type = "l", axes = FALSE,
     col = col[4], lwd = 2, xlim = xlim2, ylim = ylim3,
     xlab = "", ylab = "")

axis(4, at = seq(-0.5, 1, 0.5), col = col[4], col.axis = col[4])
text(x2, y2, ylab3, srt = -90, xpd = NA, col = col[4], cex = par()$cex.lab)

# ------------------------------------------------------------------------------
# Part 3

par(mar = c(6, 5, 1.5, 2.5))

x1 <- subsetData(filteredAnnualGBI, 2011 : 1851, "gbi")
y1 <- subsetData(filteredStackedTemperatureNGT.mid, 2011 : 1851, "stack")

x2 <- subsetData(filteredAnnualGBI, 2011 : 1871, "gbi")
y2 <- subsetData(filteredMAR, 2011 : 1871, "melt")

x3 <- subsetData(filteredStackedTemperatureNGT.mid, 2011 : 1871, "stack")
y3 <- subsetData(filteredMAR, 2011 : 1871, "melt")

# plot 1 -----
plot(x1, y1, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = c(-0.5, 1), ylim = c(-1.5, 3))

axis(1)
axis(2)
mtext(ylab3, side = 1, line = 3.25, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)
mtext("e", side = 3, adj = 0.05, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab * par()$cex)

points(x1, y1, pch = 19)

# plot 2 -----
plot(x2, y2, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = c(-0.5, 1), ylim = c(-100, 150))

axis(1)
axis(2)
mtext(ylab3, side = 1, line = 3.25, cex = par()$cex.lab * par()$cex)
mtext(ylab2, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)
mtext("f", side = 3, adj = 0.05, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab * par()$cex)

points(x2, y2, pch = 19)

# plot 3 -----
plot(x3, y3, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = c(-1.5, 3), ylim = c(-100, 150))

axis(1)
axis(2)
mtext(ylab1, side = 1, line = 3.25, cex = par()$cex.lab * par()$cex)
mtext(ylab2, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)
mtext("g", side = 3, adj = 0.05, padj = 0.5,
      line = -1, font = 2, cex = par()$cex.lab * par()$cex)

points(x3, y3, pch = 19)
lines(x <- seq(-1, 3, 0.1), coef(lm(y3 ~ x3 + 0)) * x,
      col = "dodgerblue4", lwd = 2)

# ------------------------------------------------------------------------------
dev.off()
