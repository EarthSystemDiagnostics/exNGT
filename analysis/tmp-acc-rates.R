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
                     name = "B18_12", delim = "\t") {

  readr::read_delim(file.path(path, file), col_types = readr::cols(),
                    delim = delim) %>%
    setNames(c("Year", name))

}

dev.compile <- function() {

  list(b18.12 = dev.read(),
       b21.12 = dev.read(file = "B21_2012_AccmRate.txt", name = "B21_12"),
       b23.12 = dev.read(file = "B23_2012_AccmRate.txt", name = "B23_12"),
       b26.12 = dev.read(file = "B26_2011_AccmRate.txt", name = "B26_12"),
       ngrip.12 = dev.read(file = "NGRIP_2012_AccmRate.txt", name = "NGRIP_12",
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
                      file = "B29_Schwager_Accmrate_we.txt", name = "B29"),
       neem = dev.neem()
       ) %>%
    purrr::reduce(dplyr::full_join, by = "Year") %>%
    dplyr::filter(Year <= 2011) %>%
    dplyr::mutate(dplyr::across(!Year, ~ .x * 1.e3)) %>% # in mm w.eq.
    as.data.frame()

}

dev.neem <- function(path = "data-raw/in-other", file = "NEEM_Acc.xlsx",
                     name = "NEEM") {

  require(readxl)

  tmp <- list(readxl::read_xlsx(path = file.path(path, file), range = "A2:B157") %>%
                setNames(c("Year", "NEEM2010S2")),
              readxl::read_xlsx(path = file.path(path, file), range = "D2:E284") %>%
                setNames(c("Year", "NEEM2008S3")),
              readxl::read_xlsx(path = file.path(path, file), range = "G2:H267") %>%
                setNames(c("Year", "NEEM2007S3")),
              readxl::read_xlsx(path = file.path(path, file), range = "J2:K156") %>%
                setNames(c("Year", "NEEM2008S2"))) %>%
    purrr::reduce(dplyr::full_join, by = "Year") %>%
    as.data.frame()

  data.frame(Year = tmp$Year, NEEM = rowMeans(tmp[, -1], na.rm = TRUE))

}

processAccumulationNGT <- function() {

  dev.compile() %>%
    makeAnomalies()
}


# ------------------------------------------------------------------------------
# Create and plot stack

filter.window <- 11

filteredStackedNGTacc <- processAccumulationNGT() %>%
  filterData(window = filter.window) %>%
  stackAllCores() %>%
  dplyr::filter(Year >= 1500)

stackedNGT <- processNGT() %>%
  stackNGT()

filteredStackedNGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()

xlab  <- "Year CE"
ylab1 <- "Accumulation anomaly (mm w.eq.)"
ylab2 <- grfxtools::LabelAxis(suffix = "anomaly")
cols <- c("darkblue", "black")
fsc <- 0.95

grfxtools::Quartz(file = "./zzz/ngt_acc_stack_vs_iso_stack.pdf")
grfxtools::Quartz(height = 6.5)
par(mar = c(5, 5, 2.5, 5))

plot(filteredStackedNGTacc, type = "l", axes = FALSE, lwd = 2, col = cols[1],
     xlim = c(1500, 2020), ylim = c(-50, 30), xlab = "", ylab = "")

axis(1)
axis(2, at = seq(-10, 30, 10), , col = cols[1], col.axis = cols[1])

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3.25, col = cols[1],
      at = 10, cex = fsc * par()$cex.lab * par()$cex, las = 0)

par(new = TRUE)

plot(filteredStackedNGT, type = "l", axes = FALSE, lwd = 2, col = cols[2],
     xlim = c(1500, 2020), ylim = c(-1.5, 3), xlab = "", ylab = "")

axis(side = 4, at = -1 : 2)
text(2120, 0.5, ylab2, srt = -90, xpd = NA,
     cex = fsc * par()$cex.lab * par()$cex, col = cols[2])

dev.off()

# ------------------------------------------------------------------------------
# Scatter plots differentiating between time periods

# subset time periods
t0 <- 2011 : 1500
t1 <- 1800 : 1500
t2 <- 1960 : 1801
t3 <- 2011 : 1961

x0 <- subsetData(filteredStackedNGT, t0, "stack")
y0 <- subsetData(filteredStackedNGTacc, t0, "stack")

x1 <- subsetData(filteredStackedNGT, t1, "stack")
y1 <- subsetData(filteredStackedNGTacc, t1, "stack")

x2 <- subsetData(filteredStackedNGT, t2, "stack")
y2 <- subsetData(filteredStackedNGTacc, t2, "stack")

x3 <- subsetData(filteredStackedNGT, t3, "stack")
y3 <- subsetData(filteredStackedNGTacc, t3, "stack")

# scatter
grfxtools::Quartz()#file = "./zzz/ngt_acc_iso_scatter_1500-2011.pdf")

plot(x1, y1, type = "n", xlab = ylab2, ylab = ylab1,
     xlim = c(-1, 2), ylim = c(-10, 30))

points(x1, y1, pch = 19, col = "black")
points(x2, y2, pch = 19, col = "dodgerblue4")
points(x3, y3, pch = 19, col = "firebrick4")

dev.off()

# overall:
cor.test(x0, y0) # 0.23 (p << 0.01)

cor.test(x1, y1) # 0.04 (p ~ 0.5)
cor.test(x2, y2) # 0.43 (p << 0.01)
cor.test(x3, y3) # 0.78 (p << 0.01)

# overall
res <- estimateCorrelation(stackedNGT, filteredStackedNGT,
                           filteredStackedNGTacc, filter.window = filter.window,
                           analysis.period = t0, nmc = 1000)
sprintf("r = %1.2f (p = %1.4f)", res$r, res$p)
# 1500--1800
res <- estimateCorrelation(stackedNGT, filteredStackedNGT,
                           filteredStackedNGTacc, filter.window = filter.window,
                           analysis.period = t1, nmc = 1000)
sprintf("r = %1.2f (p = %1.4f)", res$r, res$p)
#1801--1960
res <- estimateCorrelation(stackedNGT, filteredStackedNGT,
                           filteredStackedNGTacc, filter.window = filter.window,
                           analysis.period = t2, nmc = 1000)
sprintf("r = %1.2f (p = %1.4f)", res$r, res$p)
#1961-2011
res <- estimateCorrelation(stackedNGT, filteredStackedNGT,
                           filteredStackedNGTacc, filter.window = filter.window,
                           analysis.period = t3, nmc = 1000)
sprintf("r = %1.2f (p = %1.4f)", res$r, res$p)

# ------------------------------------------------------------------------------
# plot histogram

plotHistogram(piPeriod = 1500 : 1800, type = "anomaly", data.source = "acc",
              stack.method = "stack_all", plot.legend = FALSE,
              breaks = seq(-20, 40, 2.5), ylim = c(0, 0.1))

# ==============================================================================
# preliminary codes

# ------------------------------------------------------------------------------
# plot singe records

filter.window <- 11

ngt <- dev.compile() %>%
  makeAnomalies()

filteredNGT <- ngt %>%
  filterData(window = filter.window)

yoff <- 0.12
grfxtools::Quartz(height = 12)#, file = "./zzz/ngt_acc_records.pdf")

plot(ngt$Year, ngt$`B18_12`, type = "n", ylim = c(-0.05, 1.5),
     xlab = "Year CE", ylab = "Accumulation anomaly (a.u.)")

lines(ngt$Year, ngt$`B18_12`, lwd = 1, col = 1)
lines(ngt$Year, filteredNGT$`B18_12`, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B21_12` + 1 * yoff, lwd = 1, col = 2)
lines(ngt$Year, filteredNGT$`B21_12` + 1 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B23_12` + 2 * yoff, lwd = 1, col = 3)
lines(ngt$Year, filteredNGT$`B23_12` + 2 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B26_12` + 3 * yoff, lwd = 1, col = 4)
lines(ngt$Year, filteredNGT$`B26_12` + 3 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`NGRIP_12` + 4 * yoff, lwd = 1, col = 5)
lines(ngt$Year, filteredNGT$`NGRIP_12` + 4 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`NEEM` + 5 * yoff, lwd = 1, col = 6)
lines(ngt$Year, filteredNGT$`NEEM` + 5 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B16` + 6 * yoff, lwd = 1, col = 7)
lines(ngt$Year, filteredNGT$`B16` + 6 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B18` + 7 * yoff, lwd = 1, col = 8)
lines(ngt$Year, filteredNGT$`B18` + 7 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B21` + 8 * yoff, lwd = 1, col = 9)
lines(ngt$Year, filteredNGT$`B21` + 8 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B26` + 9 * yoff, lwd = 1, col = 10)
lines(ngt$Year, filteredNGT$`B26` + 9 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B29` + 10 * yoff, lwd = 1, col = 11)
lines(ngt$Year, filteredNGT$`B29` + 10 * yoff, lwd = 2, col = "darkgrey")

legend("bottomleft", col = 1 : 11, lwd = 1.5, lty = 1, bty = "n",
       legend = c("B18-2012", "B21-2012", "B23-2012", "B26-2012", "NGRIP-2012",
                  "NEEM", "B16", "B18", "B21", "B26", "B29"))

dev.off()

# ------------------------------------------------------------------------------
# plot pairs

yoff <- 0.2

grfxtools::Quartz(height = 9)

plot(ngt$Year, ngt$`B18_12`, type = "n", ylim = c(-0.05, 1.25),
     xlab = "Year CE", ylab = "Accumulation anomaly (a.u.)")

lines(ngt$Year, ngt$`B18` + 1 * yoff, lwd = 1, col = 1)
lines(ngt$Year, filteredNGT$`B18` + 1 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B18_12` + 2 * yoff, lwd = 1, col = 1)
lines(ngt$Year, filteredNGT$`B18_12` + 2 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B21` + 3 * yoff, lwd = 1, col = 2)
lines(ngt$Year, filteredNGT$`B21` + 3 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B21_12` + 4 * yoff, lwd = 1, col = 2)
lines(ngt$Year, filteredNGT$`B21_12` + 4 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B26` + 5 * yoff, lwd = 1, col = 4)
lines(ngt$Year, filteredNGT$`B26` + 5 * yoff, lwd = 2, col = "darkgrey")

lines(ngt$Year, ngt$`B26_12` + 6 * yoff, lwd = 1, col = 4)
lines(ngt$Year, filteredNGT$`B26_12` + 6 * yoff, lwd = 2, col = "darkgrey")

legend("bottomleft", col = c(1, 2, 4), lwd = 1.5, lty = 1, bty = "n",
       legend = c("B18", "B21", "B26"))
