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

# ------------------------------------------------------------------------------
# Figure 01 - time series and map

Quartz(file = "./fig/main-figure01-ac.pdf")

makeFigure01(panel = "ts")
dev.off()

Quartz(file = "./fig/main-figure01-b-raw.pdf", height = 6, width = 6)

makeFigure01(panel = "map")
dev.off()

# ------------------------------------------------------------------------------
# Figure 02 - histograms

Quartz(file = "./fig/main-figure02.pdf", height = 5, width = 12.5)
layout(matrix(1 : 3, 1, 3), widths = c(0.41, 0.41, 0.18))
par(cex = 1)

plotHistogram(type = "anomaly", plot.legend = FALSE)
mtext("a", side = 3, line = -0.5, cex = par()$cex.lab * par()$cex,
      las = 1, font = 2, adj = 0.02, padj = 0.5)

plotHistogram(type = "trend", breaks = seq(-0.2, 0.2, 0.02), ylim = c(0, 10),
              plot.legend = FALSE)
mtext("b", side = 3, line = -0.5, cex = par()$cex.lab * par()$cex,
      las = 1, font = 2, adj = 0.02, padj = 0.5)

par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

lg <- c("Pre-ind. distribution\n(1000-1800 CE)", "Gaussian fit",
        "2.5, 50, 97.5 %\nquantiles", "Recent value\n(2001-2011 CE)")

legend("topleft", legend = lg, lty = c(1, 1, 5, 1), lwd = c(10, 2, 1, 2.5),
       col = c(adjustcolor("black", 0.2), "black", "black", "firebrick4"),
       bty = "n", adj = c(0, 0.75), y.intersp = 1.5, seg.len = 2.5)

dev.off()

# ------------------------------------------------------------------------------
# Figure 03 - temperature time series and spectra

Quartz(file = "./fig/main-figure03.pdf", height = 5, width = 16)

makeFigure03()
dev.off()
