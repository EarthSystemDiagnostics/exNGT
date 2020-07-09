##
## aim:
## compare the NGT stack with Greenlandic instrumental temperature records.
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

# load NGT
NGT <- processNGT()
# load DMI and calculate anomalies
DMI <- readDMI() %>%
  makeAnomalies()

# ------------------------------------------------------------------------------
# Stack and filter

filter.window <- 11

# TODO ! temporary; replace by final stack version used in main part !
NGT.stack <- NGT %>%
  stackAllCores() %>%
  filterData(window = filter.window)

DMI.filter <- DMI %>%
  filterData(window = filter.window)

# ------------------------------------------------------------------------------
# Plot

col <- c("black", "#1b9e77", "#d95f02", "#7570b3")

xlab  <- "Year CE"
ylab1 <- bquote(delta^{"18"} * "O anomaly (\u2030)")
ylab2 <- bquote("Temperature anomaly" * " (" * degree * "C)")

xlim <- c(1900, 2020)
ylim1 <- c(-1, 1.5)
ylim2 <- c(-2, 3)

x <- 2040
y <- mean(ylim2)

Quartz(file = "./fig/ngt-dmi-comparison.pdf")
par(mar = c(5, 5, 0.5, 5))

plot(NGT.stack, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim1)

axis(1)
axis(2)

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

lines(NGT.stack, col = col[1], lwd = 2.5)

par(new = TRUE)

plot(NGT.stack, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim2)

axis(4)

text(x, y, ylab2, srt = -90, xpd = NA, cex = par()$cex.lab * par()$cex)

for (i in 2 : 4) lines(DMI.filter$Year, DMI.filter[, i], col = col[i], lwd = 2.5)

legend("topleft", c("NGT stack", colnames(DMI)[-1]),
       col = col, lwd = 2.5, bty = "n")

dev.off()

# ------------------------------------------------------------------------------
# Calculate correlations

period <- 1900 : 2011

i <- match(period, NGT.stack$Year)
j <- match(period, DMI.filter$Year)

cor(NGT.stack$stack[i], DMI.filter$Pituffik[j], use = "pair")
cor(NGT.stack$stack[i], DMI.filter$Upernavik[j], use = "pair")
cor(NGT.stack$stack[i], DMI.filter$Danmarkshavn[j], use = "pair")

