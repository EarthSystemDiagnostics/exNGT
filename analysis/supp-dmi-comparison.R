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

# ------------------------------------------------------------------------------
# Load NGT and DMI data

NGT <- processNGT()
DMI <- readDMI() %>%
  makeAnomalies()

# ------------------------------------------------------------------------------
# Stack and filter

filter.window <- 11

# filter and stack NGT
filteredStackedNGT <- NGT %>%
  filterData(window = filter.window) %>%
  stackNGT()

# filter DMI records
filteredDMI <- DMI %>%
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

grfxtools::Quartz(file = "./fig/ngt-dmi-comparison.pdf",
                  height = 4.5 , width = 8.9, mar = c(5, 5, 0.5, 5))

plot(filteredStackedNGT, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim1)

axis(1)
axis(2)

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

mtext("a", side = 3, adj = 0.01, padj = -0.35,
      line = -1, font = 2, cex = par()$cex.lab)
mtext("b", side = 3, adj = 0.99, padj = -0.35,
      line = -1, font = 2, cex = par()$cex.lab)

lines(filteredStackedNGT, col = col[1], lwd = 2.5)

par(new = TRUE)

plot(filteredStackedNGT, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = xlim, ylim = ylim2)

axis(4)

text(x, y, ylab2, srt = -90, xpd = NA, cex = par()$cex.lab * par()$cex)

for (i in 2 : 4) lines(filteredDMI$Year, filteredDMI[, i], col = col[i], lwd = 2.5)

legend("topleft", c("NGT stack", colnames(DMI)[-1]),
       col = col, lwd = 2.5, bty = "n", inset = c(0, 0.02))

dev.off()

# ------------------------------------------------------------------------------
# Calculate correlations

period <- 1900 : 2011

i <- match(period, filteredStackedNGT$Year)
j <- match(period, filteredDMI$Year)

cor(filteredStackedNGT$stack[i], filteredDMI$Pituffik[j], use = "pair")
cor(filteredStackedNGT$stack[i], filteredDMI$Upernavik[j], use = "pair")
cor(filteredStackedNGT$stack[i], filteredDMI$Danmarkshavn[j], use = "pair")

