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

snrNGT <- analysisRecords %>%
  proxysnr::ArraySpectra(df.log = 0.1) %>%
  proxysnr::SeparateSpectra() %>%
  proxysnr::StackCorrelation(N = 10,
                             freq.cut.lower = 1 / 200, freq.cut.upper = 1 / 5)

yat.snr <- c(0.1, 0.2, 0.5, 1, 2)
yat.r   <- seq(0.3, 0.8, 0.1)

Quartz(file = "./fig/ngt-snr.pdf")
par(mar = c(5, 5, 0.5, 5))

PaleoSpec::LPlot(snrNGT$signal, bPeriod = TRUE, bNoPlot = TRUE, axes = FALSE,
                 xlab = "", ylab = "", ylim = c(0.1, 2))
axis(1)
axis(2, at = yat.snr)
axis(4, at = yat.r^2 / (1 - yat.r^2), labels = yat.r)

mtext("Time period (yr)", 1, 3.5, cex = par()$cex.lab)
mtext("Signal-to-Noise Ratio", 2, 3.5, cex = par()$cex.lab, las = 0)
text(2.4, sqrt(0.1 * 2), "Correlation", xpd = NA, srt = -90, cex = par()$cex.lab)

lines(1 / snrNGT$signal$freq, snrNGT$signal$spec / snrNGT$noise$spec, lwd = 2)

dev.off()

