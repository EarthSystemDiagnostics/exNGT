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
# -> 1505-1979; all except NEGIS and NEEM
year.upper <- max(recordOccurrence$Year[which(!is.na(recordOccurrence$GRIP))])
year.lower <- min(recordOccurrence$Year[which(!is.na(recordOccurrence$B26))])

# visualize
grfxtools::Quartz(mar = c(5, 5, 0.5, 1))
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

# ------------------------------------------------------------------------------
# Do the SNR analysis

# select time and record range according to above and apply spectral analysis
spectraNGT <- selectNGTForSpectra() %>%
  proxysnr::ArraySpectra(df.log = 0.1) %>%
  proxysnr::SeparateSpectra()

# min, max and mean number of records contributing to NGT stack 1000-2011 CE
N <- processNGT() %>%
  stackNGT(stack = FALSE) %>%
  countRecords() %>%
  dplyr::filter(Year >= 1000) %>%
  dplyr::summarise(min = min(n), max = max(n), mean = round(mean(n)))

# obtain signal-to-noise ratio spectrum between 1/200 and 1/5 yr^-1
f1 <- 1 / 200
f2 <- 1 / 5
snrNGT.min <- spectraNGT %>%
  proxysnr::IntegratedSNR(N = N$min, freq.cut.lower = f1, freq.cut.upper = f2)
snrNGT.max <- spectraNGT %>%
  proxysnr::IntegratedSNR(N = N$max, freq.cut.lower = f1, freq.cut.upper = f2)
snrNGT.mean <- spectraNGT %>%
  proxysnr::IntegratedSNR(N = N$mean, freq.cut.lower = f1, freq.cut.upper = f2)

# obtain timescale threshold below which spectra are dominated by diffusion;
# i.e. timescale at which diffusion transfer function drops to 1/e given the
# maximum estimated diffusion length in years across all NGT records
getDiffusionThreshold <- function(s, k) {
  if (k >= 1 | k <= 0) stop("Need 0 < k < 1.", call. = FALSE)
  sqrt(-1 * (2 * pi *s)^2 / log(k))
}
diffThreshold <- diffuseNGT(NGT, result = "sigma") %>%
  .$sigma %>%
  dplyr::select(-Time) %>%
  max() %>%
  getDiffusionThreshold(k = exp(-1))

# ------------------------------------------------------------------------------
# Plot

yat.snr <- c(0.5, 1, 2, 5, 10, 20)
yat.r   <- c(seq(0.6, 0.9, 0.1), 0.95, 0.97)

xlab  <- "Time period (yr)"
ylab1 <- grfxtools::LabelAxis("Signal PSD", time.unit = "yr", unit.type = "psd")
ylab2 <- "Signal-to-Noise Ratio"
ylab3 <- "Correlation"

xlim <- c(225, 5)
ylim1 <- c(0.05, 5)
ylim2 <- c(0.5, 20)

label1 <- expression(bold("a"))
label2 <- expression(bold("b"))

n.crit.lower <- 1
n.crit.upper <- length(which(spectraNGT$signal$freq > f2))

grfxtools::Quartz(file = "./fig/supplement-ngt-snr.pdf",
                  height = 4.5, width = 12, mfcol = c(1, 2),
                  oma = c(5, 0, 0.5, 5.5), mar = c(0, 5.5, 0, 0))

proxysnr:::LPlot(spectraNGT$signal, bPeriod = TRUE, bNoPlot = TRUE, axes = FALSE,
                 xlab = "", ylab = "", xlim = xlim, ylim = ylim1)

axis(1)
axis(2)

mtext(xlab, 1, 3.5, cex = par()$cex.lab)
mtext(ylab1, 2, 3.75, cex = par()$cex.lab, las = 0)
mtext(label1, side = 3, line = -1.5, cex = par()$cex.lab,
      adj = 0.02, padj = 0.2)

grfxtools::Polyplot(c(diffThreshold, 5), rep(ylim1[1], 2), rep(ylim1[2], 2),
                    col = "firebrick4")

proxysnr:::LLines(spectraNGT$signal, bPeriod = TRUE, lwd = 2,
                  removeFirst = n.crit.lower, removeLast = n.crit.upper)

proxysnr:::LPlot(spectraNGT$signal, bPeriod = TRUE, bNoPlot = TRUE, axes = FALSE,
                 xlab = "", ylab = "", xlim = xlim, ylim = ylim2)
axis(1)
axis(2, at = yat.snr)
axis(4, at = yat.r^2 / (1 - yat.r^2), labels = sprintf("%1.2f", yat.r))

mtext(xlab, 1, 3.5, cex = par()$cex.lab)
mtext(ylab2, 2, 3.5, cex = par()$cex.lab, las = 0)
text(1.85, sqrt(prod(range(yat.snr))), ylab3, xpd = NA, srt = -90,
     cex = par()$cex.lab)
mtext(label2, side = 3, line = -1.5, cex = par()$cex.lab,
      adj = 0.02, padj = 0.2)

grfxtools::Polyplot(c(diffThreshold, 5), rep(ylim2[1], 2), rep(ylim2[2], 2),
                    col = "firebrick4")

grfxtools::Polyplot(1 / snrNGT.min$freq, snrNGT.min$spec, snrNGT.max$spec)
proxysnr:::LLines(snrNGT.min, bPeriod = TRUE, lwd = 1)
proxysnr:::LLines(snrNGT.max, bPeriod = TRUE, lwd = 1)
proxysnr:::LLines(snrNGT.mean, bPeriod = TRUE, lwd = 2)

dev.off()
