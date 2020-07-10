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
# Process

filter.window <- 11
permil2temperature <- 1 / 0.67

# NGT: stack and convert to temperature
# TODO ! temporary; replace by final stack version !
NGT.stack <- NGT %>%
  stackAllCores() %>%
  dplyr::mutate(stack.temperature = permil2temperature * stack)

# filter NGT
NGT.stack.f <- NGT.stack %>%
  filterData(window = filter.window)

# filter Arctic2k
Arctic2k.f <- Arctic2k %>%
  filterData(window = filter.window)

# ------------------------------------------------------------------------------
# Fig. 1, main, full NGT

xlab  <- "Year CE"
ylab <- bquote(delta^{"18"} * "O anomaly (\u2030)")

xlim <- c(1000, 2020)
ylim <- c(-2.5, 2.5)

startNew <- 1993
i <- match(startNew, NGT.stack$Year)
n <- nrow(NGT.stack)

Quartz(file = "./fig/main-01-ngt.pdf", height = 4.5)

plot(NGT.stack$Year, NGT.stack$stack,
     type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim)

axis(1)
axis(2)

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

lines(NGT.stack$Year, NGT.stack$stack, col = "darkgrey")

lines(NGT.stack.f$Year[i : n], NGT.stack.f$stack[i : n],
      col = "black", lwd = 2.5)
lines(NGT.stack.f$Year[1 : i], NGT.stack.f$stack[1 : i],
      col = "firebrick3", lwd = 2.5)

dev.off()

# ------------------------------------------------------------------------------
# Fig. 1, minor, NGT and Arctic2k

col <- c("black", "dodgerblue4")

xlab  <- "Year CE"
ylab <- bquote("Temperature anomaly" * " (" * degree * "C)")

xlim <- c(1840, 2020)
ylim <- c(-3, 3)

Quartz(file = "./fig/main-02-ngt-arctic2k-comparison.pdf", height = 4.5)

plot(NGT.stack$Year, NGT.stack$stack,
     type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim)

axis(1, at = seq(xlim[1], xlim[2], 20))
axis(2)

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

lines(NGT.stack$Year, NGT.stack$stack, col = "darkgrey")

lines(NGT.stack.f$Year, NGT.stack.f$stack, col = col[1], lwd = 2.5)

lines(Arctic2k$Year, Arctic2k$TempAnomaly,
      col = adjustcolor(col[2], alpha = 0.6))
lines(Arctic2k.f$Year, Arctic2k.f$TempAnomaly, col = col[2], lwd = 2.5)

legend("topleft", c("NGT stack", "Arctic2k"), col = col, lwd = 2.5, bty = "n")

dev.off()

