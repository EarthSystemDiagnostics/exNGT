##
## aim:
## preliminary code for accumulation rate analyses.
## relation:
## NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

# ------------------------------------------------------------------------------
# compile data

dev.read <- function(path = "data-raw/in-raw", file = "B18_2012_AccmRate.txt",
                     name = "B18-2012", delim = "\t") {

  readr::read_delim(file.path(path, file), col_types = readr::cols(),
                    delim = delim) %>%
    setNames(c("Year", name))

}

dev.compile <- function() {

  list(b18.12 = dev.read(),
       b21.12 = dev.read(file = "B21_2012_AccmRate.txt", name = "B21-2012"),
       b23.12 = dev.read(file = "B23_2012_AccmRate.txt", name = "B23-2012"),
       ngrip.12 = dev.read(file = "NGRIP_2012_AccmRate.txt", name = "NGRIP-2012",
                        delim = " "),
       b16 = dev.read(path = "data-raw/in-other",
                      file = "B16_Schwager_Accmrate_we.txt", name = "B16"),
       b18 = dev.read(path = "data-raw/in-other",
                      file = "B18_Schwager_Accmrate_we.txt", name = "B18"),
       b21 = dev.read(path = "data-raw/in-other",
                      file = "B21_Schwager_Accmrate_we.txt", name = "B21"),
       b26 = dev.read(path = "data-raw/in-other",
                      file = "B26_Schwager_Accmrate_we.txt", name = "B26"),
       b29 = dev.read(path = "data-raw/in-other",
                      file = "B29_Schwager_Accmrate_we.txt", name = "B29")
       ) %>%
    purrr::reduce(dplyr::full_join, by = "Year") %>%
    dplyr::filter(Year <= 2011) %>%
    as.data.frame()
 
}

# ------------------------------------------------------------------------------
# plot singe records

filter.window <- 11

ngt <- dev.compile() %>%
  makeAnomalies()

filteredNGT <- ngt %>%
  filterData(window = filter.window)

yoff <- 0.12
grfxtools::Quartz(height = 10, file = "./zzz/ngt_acc_records.pdf")

plot(ngt$Year, ngt$`B18-2012`, type = "n", ylim = c(-0.05, 1),
     xlab = "Year CE", ylab = "Accumulation anomaly (a.u.)")

lines(ngt$Year, ngt$`B18-2012`, lwd = 1, col = 1)
lines(ngt$Year, filteredNGT$`B18-2012`, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B21-2012` + 1 * yoff, lwd = 1, col = 2)
lines(ngt$Year, filteredNGT$`B21-2012` + 1 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B23-2012` + 2 * yoff, lwd = 1, col = 3)
lines(ngt$Year, filteredNGT$`B23-2012` + 2 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`NGRIP-2012` + 3 * yoff, lwd = 1, col = 4)
lines(ngt$Year, filteredNGT$`NGRIP-2012` + 3 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B16` + 4 * yoff, lwd = 1, col = 5)
lines(ngt$Year, filteredNGT$`B16` + 4 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B18` + 5 * yoff, lwd = 1, col = 6)
lines(ngt$Year, filteredNGT$`B18` + 5 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B21` + 6 * yoff, lwd = 1, col = 7)
lines(ngt$Year, filteredNGT$`B21` + 6 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B26` + 7 * yoff, lwd = 1, col = 8)
lines(ngt$Year, filteredNGT$`B26` + 7 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B29` + 8 * yoff, lwd = 1, col = 9)
lines(ngt$Year, filteredNGT$`B29` + 8 * yoff, lwd = 2, col = "darkgrey")

legend("bottomleft", col = 1 : 9, lwd = 1.5, lty = 1, bty = "n",
       legend = c("B18-2012", "B21-2012", "B23-2012", "NGRIP-2012",
                  "B16", "B18", "B21", "B26", "B29"))

dev.off()

# ------------------------------------------------------------------------------
# create and plot stack

filteredAccStack <- data.frame(Year = ngt$Year,
                               acc = rowMeans(filteredNGT[, -1], na.rm = TRUE))

filteredStackedNGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()

grfxtools::Quartz(file = "./zzz/ngt_acc_stack_vs_iso_stack.pdf")
par(mar = c(5, 5, 0.5, 5))

plot(filteredAccStack, type = "l", lwd = 2, col = 1,
     xlim = c(1500, 2020), ylim = c(-0.04, 0.04),
     xlab = "Year CE", ylab = "Accumulation anomaly (m w.eq.)")

par(new = TRUE)

plot(filteredStackedNGT, type = "l", axes = FALSE, lwd = 2, col = 2,
     xlim = c(1500, 2020), ylim = c(-2, 2),
     xlab = "", ylab = "")
axis(side = 4, col = 2)
mtext(grfxtools::LabelAxis(suffix = "anomaly"), side = 4, col = 2,
      line = 3.5, las = 0, cex = par()$cex.lab)

dev.off()

# subset time period 1500 - 2011
t <- 2011 : 1500

x <- subsetData(filteredStackedNGT, t, "stack")
y <- subsetData(filteredAccStack, t, "acc")

# scatter
grfxtools::Quartz(file = "./zzz/ngt_acc_iso_scatter_1500-2011.pdf")
plot(x, y, pch = 19,
     xlab = grfxtools::LabelAxis(), ylab = "Accumulation (m w.eq.)")
dev.off()

cor.test(x, y) # 0.24 (p << 0.01)

# subset time period 1800 - 1500
t <- 1800 : 1500

x <- subsetData(filteredStackedNGT, t, "stack")
y <- subsetData(filteredAccStack, t, "acc")

# scatter
grfxtools::Quartz(file = "./zzz/ngt_acc_iso_scatter_1500-1800.pdf")
plot(x, y, pch = 19,
     xlab = grfxtools::LabelAxis(), ylab = "Accumulation (m w.eq.)")
dev.off()
cor.test(x, y) # 0.04 (p = 0.5)


