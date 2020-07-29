##
## aim:
## analyse timescale-dependent SNR of NGT records.
## relation:
## NGT paper supplementary; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

library(magrittr)
library(proxysnr)
library(PaleoSpec)

# ------------------------------------------------------------------------------
# Input

# all NGT anomaly records
NGT <- processNGT()

# NGT records of main stack at annual resolution
mainStackRecords <- NGT %>%
  stackNGT(stack = FALSE)

# ------------------------------------------------------------------------------
# Analyse record occurrences

recordOccurrence <- data.frame(
  Year = NGT$Year,
  matrix(1 : (ncol(mainStackRecords) - 1),
         nrow = nrow(NGT), ncol = ncol(mainStackRecords) - 1,
         byrow = TRUE)) %>%
  setNames(colnames(mainStackRecords))

recordOccurrence[which(is.na(mainStackRecords), arr.ind = TRUE)] <- NA

# set upper and lower years for selecting 14 records
# -> 1505-1978; all except NEGIS and NEEM
year.upper <- max(recordOccurrence$Year[which(!is.na(recordOccurrence$GRIP))])
year.lower <- min(recordOccurrence$Year[which(!is.na(recordOccurrence$B26))])

# visualize
Quartz()
par(mar = c(5, 5, 0.5, 1))
plot.new()
plot.window(xlim = c(1000, 2011), ylim = c(1, 16), xaxs = 'i')
box()

axis(1)
axis(2, at = 1 : 16, labels = colnames(recordOccurrence)[-1])
mtext("Year CE", 1, 3.5, cex = par()$cex.lab)

for (i in 2 : 17) {

  if (colnames(recordOccurrence)[i] %in% c("NEGIS", "NEEM")) {
    col <- "darkgrey"
  } else {
    col <- "black"
  }
  lines(recordOccurrence$Year, recordOccurrence[, i], lwd = 1.5, col = col)
}

abline(v = c(year.upper, year.lower), col = "darkgrey", lty = 2)

# select time and record range
i <- match(year.upper : year.lower, mainStackRecords$Year)
j <- colnames(mainStackRecords)[-match(c("Year", "NEEM", "NEGIS"),
                                       colnames(mainStackRecords))]

analysisRecords <- mainStackRecords[i, j]

# ------------------------------------------------------------------------------
# Do the SNR analysis

spectraNGT <- analysisRecords %>%
  proxysnr::ArraySpectra(df.log = 0.1) %>%
  proxysnr::SeparateSpectra()

snrNGT <- spectraNGT %>%
  proxysnr::StackCorrelation(N = 10,
                             freq.cut.lower = 1 / 200, freq.cut.upper = 1 / 5)

yat.snr <- c(0.1, 0.2, 0.5, 1, 2)
yat.r   <- seq(0.3, 0.8, 0.1)

xlab  <- "Time period (yr)"
ylab1 <- expression("Signal PSD " * "(\u2030"^{2}%.%"yr)")
ylab2 <- "Signal-to-Noise Ratio"
ylab3 <- "Correlation"

label1 <- expression(bold("a"))
label2 <- expression(bold("b"))

n.crit.lower <- 1
n.crit.upper <- length(which(spectraNGT$signal$freq > 1 / 5))

Quartz(file = "./fig/ngt-snr.pdf", height = 4.5, width = 12)
op <- par(LoadGraphicsPar(mfcol = c(1, 2), oma = c(5, 0, 0.5, 5.5),
                          mar = c(0, 5.5, 0, 0)))

PaleoSpec::LPlot(spectraNGT$signal, bPeriod = TRUE, bNoPlot = TRUE, axes = FALSE,
                 xlab = "", ylab = "", xlim = c(225, 5), ylim = c(0.05, 5))

axis(1)
axis(2)

mtext(xlab, 1, 3.5, cex = par()$cex.lab)
mtext(ylab1, 2, 3.75, cex = par()$cex.lab, las = 0)
mtext(label1, side = 3, line = -1.5, cex = par()$cex.lab,
      adj = 0.02, padj = 0.2)

PaleoSpec::LLines(spectraNGT$signal, bPeriod = TRUE, lwd = 2,
                  removeFirst = n.crit.lower, removeLast = n.crit.upper)

PaleoSpec::LPlot(snrNGT$signal, bPeriod = TRUE, bNoPlot = TRUE, axes = FALSE,
                 xlab = "", ylab = "", xlim = c(225, 5), ylim = c(0.1, 2))
axis(1)
axis(2, at = yat.snr)
axis(4, at = yat.r^2 / (1 - yat.r^2), labels = yat.r)

mtext(xlab, 1, 3.5, cex = par()$cex.lab)
mtext(ylab2, 2, 3.5, cex = par()$cex.lab, las = 0)
text(2, sqrt(0.1 * 2), ylab3, xpd = NA, srt = -90, cex = par()$cex.lab)
mtext(label2, side = 3, line = -1.5, cex = par()$cex.lab,
      adj = 0.02, padj = 0.2)

lines(1 / snrNGT$signal$freq, snrNGT$signal$spec / snrNGT$noise$spec, lwd = 2)

par(op)
dev.off()

