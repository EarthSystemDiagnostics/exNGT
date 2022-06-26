##
## aim:
## supplementary accumulation rate analyses
## relation:
## NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

# ------------------------------------------------------------------------------
# Get data and estimate overall correlation

filter.window <- 11

filteredStackedAccumulation <- processAccumulation() %>%
  filterData(window = filter.window) %>%
  stackAllCores() %>%
  dplyr::filter(Year >= 1500)

stackedNGT <- processNGT() %>%
  stackNGT() %>%
  dplyr::filter(Year >= 1500)

filteredStackedNGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT() %>%
  dplyr::filter(Year >= 1500)

# ------------------------------------------------------------------------------
# Prepare plotting of isotope-to-accumulation stack comparison

xlab  <- "Year CE"
ylab1 <- grfxtools::LabelAxis("Accumulation rate", unit = "mm w.eq.",
                              time.unit = "yr", unit.type = "trend")
ylab2 <- grfxtools::LabelAxis("NGT-2012")
cols <- c("deepskyblue4", "black")

asp <- 1
w <- 13.6
natfig(file = "./fig/supplement-ngt-2012-accumulation.pdf",
       height = asp * w, width = w, mar = c(5, 5, 5, 5))
layout(matrix(c(1, 2, 1, 3), 2, 2))
par(cex = 1)

# ------------------------------------------------------------------------------
# I. Plot stacks

plot(filteredStackedAccumulation, type = "l", axes = FALSE, lwd = Lwd(2),
     col = cols[1], xlim = c(1500, 2020), ylim = c(-50, 30),
     xlab = "", ylab = "")

axisLwd(1)
axisLwd(2, at = seq(-10, 30, 10), col = cols[1], col.axis = cols[1])

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3, col = cols[1],
      at = 10, cex = par()$cex.lab * par()$cex, las = 0)

mtext("a", side = 3, adj = 0.01, padj = 0.5, line = 2, font = 2, cex = 4 / 3)

par(new = TRUE)

plot(filteredStackedNGT, type = "l", axes = FALSE, lwd = Lwd(2), col = cols[2],
     xlim = c(1500, 2020), ylim = c(-2.5, 5), xlab = "", ylab = "")

axisLwd(side = 4, at = -2 : 2)
text(2085, 0, ylab2, srt = -90, xpd = NA,
     cex = par()$cex.lab, col = cols[2])

# ------------------------------------------------------------------------------
# II. Plot scatter plot

plot(filteredStackedNGT$stack, filteredStackedAccumulation$stack,
     type = "p", pch = 19, axes = FALSE, xlab = "", ylab = "",
     xlim = c(-2, 2), ylim = c(-10, 30))

axisLwd(1)
axisLwd(2)
mtext(ylab2, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3, cex = par()$cex.lab * par()$cex, las = 0)

mtext("b", side = 3, adj = 0.01, padj = 0.5, line = 2, font = 2, cex = 4 / 3)

# ------------------------------------------------------------------------------
# III. Plot histogram

plotHistogram(piPeriod = 1500 : 1800, type = "anomaly", data.source = "acc",
              stack.method = "stack_all", plot.legend = FALSE,
              breaks = seq(-20, 40, 2.5), ylim = c(0, 0.1))
mtext("c", side = 3, adj = 0.01, padj = 0.5, line = 2, font = 2, cex = 4 / 3)

dev.off()


# ==============================================================================
# further analyses only for review
# ==============================================================================


# ------------------------------------------------------------------------------
# plot single records

filter.window <- 11

ngt <- processAccumulation() %>%
  dplyr::filter(Year >= 1500)

filteredNGT <- ngt %>%
  filterData(window = filter.window)

yoff <- 0.14 * 1.e3
yoff2 <- 0.075 * 1.e3
cols <- c("black", "red", "blue")

grfxtools::Quartz(height = 12, file = "./zzz/ngt_all_acc_records.pdf")

plot(ngt$Year, ngt$`B18_12`, type = "n", ylim = c(-0.05, 1.5) * 1.e3,
     xlab = "Year CE", ylab = "Accumulation anomaly (a.u.)")

lines(ngt$Year, ngt$`B18_12`, lwd = 1, col = cols[1])
lines(ngt$Year, filteredNGT$`B18_12`, lwd = 3, col = "darkgrey")
text(x = 1525, y = 0, "B18-12", col = cols[1], cex = 1.25)

lines(ngt$Year, ngt$`B21_12` + 1 * yoff, lwd = 1, col = cols[2])
lines(ngt$Year, filteredNGT$`B21_12` + 1 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 1 * yoff, "B21-12", col = cols[2], cex = 1.25)

lines(ngt$Year, ngt$`B23_12` + 2 * yoff, lwd = 1, col = cols[3])
lines(ngt$Year, filteredNGT$`B23_12` + 2 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 2 * yoff, "B23-12", col = cols[3], cex = 1.25)

lines(ngt$Year, ngt$`B26_12` + 3 * yoff, lwd = 1, col = cols[1])
lines(ngt$Year, filteredNGT$`B26_12` + 3 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 3 * yoff, "B26-12", col = cols[1], cex = 1.25)

lines(ngt$Year, ngt$`NGRIP_12` + 4 * yoff, lwd = 1, col = cols[2])
lines(ngt$Year, filteredNGT$`NGRIP_12` + 4 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 4 * yoff, "NGRIP-12", col = cols[2], cex = 1.25)

lines(ngt$Year, ngt$`NEEM` + 5 * yoff, lwd = 1, col = cols[3])
lines(ngt$Year, filteredNGT$`NEEM` + 5 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 5 * yoff, "NEEM", col = cols[3], cex = 1.25)

lines(ngt$Year, ngt$`B16` + 6 * yoff, lwd = 1, col = cols[1])
lines(ngt$Year, filteredNGT$`B16` + 6 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 6 * yoff + yoff2, "B16", col = cols[1], cex = 1.25)

lines(ngt$Year, ngt$`B18` + 7 * yoff, lwd = 1, col = cols[2])
lines(ngt$Year, filteredNGT$`B18` + 7 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 7 * yoff + yoff2, "B18", col = cols[2], cex = 1.25)

lines(ngt$Year, ngt$`B21` + 8 * yoff, lwd = 1, col = cols[3])
lines(ngt$Year, filteredNGT$`B21` + 8 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 8 * yoff + yoff2, "B21", col = cols[3], cex = 1.25)

lines(ngt$Year, ngt$`B26` + 9 * yoff, lwd = 1, col = cols[1])
lines(ngt$Year, filteredNGT$`B26` + 9 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 9 * yoff + yoff2, "B26", col = cols[1], cex = 1.25)

lines(ngt$Year, ngt$`B29` + 10 * yoff, lwd = 1, col = cols[2])
lines(ngt$Year, filteredNGT$`B29` + 10 * yoff, lwd = 3, col = "darkgrey")
text(x = 1525, y = 10 * yoff + yoff2, "B29", col = cols[2], cex = 1.25)

dev.off()


# ------------------------------------------------------------------------------
# spectrum/SNR analysis similar to d18O data

filteredAccumulation <- processAccumulation() %>%
  filterData(window = filter.window)

recordOccurrence <- data.frame(
  Year = filteredAccumulation$Year,
  matrix(1 : (ncol(filteredAccumulation) - 1),
         nrow = nrow(filteredAccumulation), ncol = ncol(filteredAccumulation) - 1,
         byrow = TRUE)) %>%
  setNames(colnames(filteredAccumulation))

recordOccurrence[which(is.na(filteredAccumulation), arr.ind = TRUE)] <- NA

grfxtools::Quartz(mar = c(5, 6.5, 0.5, 1))
plot.new()
plot.window(xlim = c(1000, 2011), ylim = c(1, 11), xaxs = 'i')
box()

axis(1)
axis(2, at = 1 : 11, labels = colnames(recordOccurrence)[-1])
mtext("Year CE", 1, 3.5, cex = par()$cex.lab)

for (i in 2 : 12) {
  lines(recordOccurrence$Year, recordOccurrence[, i], lwd = 1.5, col = "black")
}

noMissingVal <- function(x) {!any(is.na(x))}
spectraNGT <- filteredAccumulation %>%
  dplyr::slice(match(1992 : 1502, Year)) %>%
  dplyr::select(where(noMissingVal)) %>%
  dplyr::select(-Year) %>%
  proxysnr::ArraySpectra(df.log = 0.1) %>%
  proxysnr::SeparateSpectra()

# min, max and mean number of records contributing to stack
N <- filteredAccumulation %>%
  countRecords() %>%
  dplyr::filter(Year >= 1502) %>%
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

yat.snr <- c(0.5, 1, 2, 5, 10, 20)

xlab  <- "Time period (yr)"
ylab1 <- grfxtools::LabelAxis("Signal PSD", time.unit = "yr", unit.type = "psd")
ylab2 <- "Signal-to-Noise Ratio"

xlim <- c(225, 5)
ylim <- c(0.5, 20)

grfxtools::Quartz(height = 4.5, width = 6, mar = c(5, 5.5, 0.5, 0.5),
                  file = "./zzz/ngt_accumulation_snr.pdf")

proxysnr:::LPlot(snrNGT.mean, bPeriod = TRUE, bNoPlot = TRUE, axes = FALSE,
                 xlab = "", ylab = "", xlim = xlim, ylim = ylim)

axis(1)
axis(2, at = yat.snr)
mtext(xlab, 1, 3.5, cex = par()$cex.lab)
mtext(ylab2, 2, 3.5, cex = par()$cex.lab, las = 0)

grfxtools::Polyplot(1 / snrNGT.min$freq, snrNGT.min$spec, snrNGT.max$spec)
proxysnr:::LLines(snrNGT.min, bPeriod = TRUE, lwd = 1)
proxysnr:::LLines(snrNGT.max, bPeriod = TRUE, lwd = 1)
proxysnr:::LLines(snrNGT.mean, bPeriod = TRUE, lwd = 2)

dev.off()

# ------------------------------------------------------------------------------
# scatter plot with non-anomaly data

path <- "data/NGT2012_Accumulation_20220520.csv"
ngt <- list(oxy = readNGT(), acc = read.csv(path, header = TRUE)) %>%
  lapply(filterData, window = 11) %>%
  lapply(function(x) {dplyr::filter(x, Year >= 1500)})

grfxtools::Quartz(height = 6, width = 6,
                  file = "./zzz/ngt_d18O_accumulation_nonanom_buchardt.pdf")

plot(ngt$oxy$B18, ngt$acc$B18, type = "n", xlim = c(-40, -30), ylim = c(0, 250),
     xlab = grfxtools::LabelAxis(), ylab = "Accumulation rate (mm w.eq.)")

for (nm in names(ngt$acc)[-1]) {
  points(ngt$oxy[, nm], ngt$acc[, nm], pch = 19)
}

# Buchardt fits
x <- seq(-40, -30, 0.1)
y.ne <- 230 * 917 / 1000 * exp(0.023 * (x + 16.75))
y.nw <- 230 * 917 / 1000 * exp(0.1 * (x + 33.64))

lines(x, y.ne, lwd = 2.5, col = "green")
lines(x, y.nw, lwd = 2.5, col = "darkgrey")

legend("bottomright", c("Buchardt-NE", "Buchardt-NW"),
       lty = 1, lwd = 2.5, col = c("green", "black"), bty = "n")

dev.off()
