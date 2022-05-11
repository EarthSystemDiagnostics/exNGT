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

  #require(readxl)

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
# Prepare plotting

filter.window <- 11

filteredAccumulationNGT <- processAccumulationNGT() %>%
  filterData(window = filter.window)

filteredStackedAccumulationNGT <- filteredAccumulationNGT %>%
  stackAllCores() %>%
  dplyr::filter(Year >= 1500)

# use the same records for stacking d18O as for accumulation
stackedNGT <- processNGT() %>%
  dplyr::select(names(filteredAccumulationNGT)) %>%
  stackAllCores() %>%
  dplyr::filter(Year >= 1500)

filteredStackedNGT <- processNGT() %>%
  dplyr::select(names(filteredAccumulationNGT)) %>%
  filterData(window = filter.window) %>%
  stackAllCores() %>%
  dplyr::filter(Year >= 1500)

# main NGT-2012 d18O stack
NGT2012 <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()

xlab  <- "Year CE"
ylab1 <- grfxtools::LabelAxis("Accumulation rate", unit = "mm w.eq.",
                              time.unit = "yr", unit.type = "trend")
ylab2 <- grfxtools::LabelAxis("NGT-2012")
ylab3 <- grfxtools::LabelAxis()
cols <- c("deepskyblue4", "black")

grfxtools::Quartz(height = 12, width = 12,
                  mar = c(5, 5, 5, 5),
                  file = "./fig/supplement-ngt-2012-accumulation.pdf")
layout(matrix(c(1, 2, 1, 3), 2, 2))

# ------------------------------------------------------------------------------
# I. Plot stack

plot(filteredStackedAccumulationNGT, type = "l", axes = FALSE, lwd = 2,
     col = cols[1], xlim = c(1500, 2020), ylim = c(-50, 30),
     xlab = "", ylab = "")

axis(1)
axis(2, at = seq(-10, 30, 10), col = cols[1], col.axis = cols[1])

mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3.25, col = cols[1],
      at = 10, cex = par()$cex.lab * par()$cex, las = 0)

mtext("a", side = 3, adj = 0.01, padj = 0.5,
      line = 2, font = 2, cex = par()$cex.lab)

par(new = TRUE)

plot(NGT2012, type = "l", axes = FALSE, lwd = 2, col = cols[2],
     xlim = c(1500, 2020), ylim = c(-2.5, 5), xlab = "", ylab = "")

axis(side = 4, at = -2 : 2)
text(2075, 0, ylab2, srt = -90, xpd = NA,
     cex = par()$cex.lab, col = cols[2])

# ------------------------------------------------------------------------------
# II. Scatter plots differentiating between time periods

# subset time periods
t0 <- 2011 : 1500
t1 <- 1800 : 1500
t2 <- 1960 : 1801
t3 <- 2011 : 1961

x0 <- subsetData(filteredStackedNGT, t0, "stack")
y0 <- subsetData(filteredStackedAccumulationNGT, t0, "stack")

x1 <- subsetData(filteredStackedNGT, t1, "stack")
y1 <- subsetData(filteredStackedAccumulationNGT, t1, "stack")

x2 <- subsetData(filteredStackedNGT, t2, "stack")
y2 <- subsetData(filteredStackedAccumulationNGT, t2, "stack")

x3 <- subsetData(filteredStackedNGT, t3, "stack")
y3 <- subsetData(filteredStackedAccumulationNGT, t3, "stack")

# get correlations
# overall
res0 <- estimateCorrelation(stackedNGT, filteredStackedNGT,
                            filteredStackedAccumulationNGT,
                            filter.window = filter.window,
                            analysis.period = t0, nmc = 10000)
sprintf("r = %1.2f (p = %1.2f)", res0$r, res0$p) # 0.23 (p = 0.06)
# 1500--1800
res1 <- estimateCorrelation(stackedNGT, filteredStackedNGT,
                            filteredStackedAccumulationNGT,
                            filter.window = filter.window,
                            analysis.period = t1, nmc = 10000)
sprintf("r = %1.2f (p = %1.2f)", res1$r, res1$p) # 0.02 (p = 0.47)
#1801--1960
res2 <- estimateCorrelation(stackedNGT, filteredStackedNGT,
                            filteredStackedAccumulationNGT,
                            filter.window = filter.window,
                            analysis.period = t2, nmc = 10000)
sprintf("r = %1.2f (p = %1.2f)", res2$r, res2$p) # 0.44 (p = 0.05)
#1961-2011
res3 <- estimateCorrelation(stackedNGT, filteredStackedNGT,
                            filteredStackedAccumulationNGT,
                            filter.window = filter.window,
                            analysis.period = t3, nmc = 10000)
sprintf("r = %1.2f (p = %1.2f)", res3$r, res3$p) # 0.78 (p = 0.06)

# scatter plot
plot(x1, y1, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(-2, 2), ylim = c(-10, 30))

axis(1)
axis(2)
mtext(ylab3, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(ylab1, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

mtext("b", side = 3, adj = 0.01, padj = 0.5,
      line = 2, font = 2, cex = par()$cex.lab)

points(x1, y1, pch = 19, col = "black")
points(x2, y2, pch = 19, col = "black")
points(x3, y3, pch = 19, col = "black")

## legend("bottomright",
##        c(sprintf("1500-1800 CE (%1.2f, p = %1.2f)", res1$r, res1$p),
##          sprintf("1801-1960 CE (%1.2f, p = %1.2f)", res2$r, res2$p),
##          sprintf("1961-2011 CE (%1.2f, p = %1.2f)", res3$r, res3$p)),
##        pch = 19, col = c("black", "dodgerblue4", "firebrick4"),
##        inset = c(0.025, 0), bty = "n")

# ------------------------------------------------------------------------------
# III. Plot histogram

plotHistogram(piPeriod = 1500 : 1800, type = "anomaly", data.source = "acc",
              stack.method = "stack_all", plot.legend = FALSE,
              breaks = seq(-20, 40, 2.5), ylim = c(0, 0.1))
mtext("c", side = 3, adj = 0.01, padj = 0.5,
      line = 2, font = 2, cex = par()$cex.lab)
dev.off()

# ==============================================================================
# preliminary codes

selectFixedNumber <- function(NGT, nfix) {

  rslt <- NULL
  nrow <- nrow(NGT)

  for (i in 1 : nrow) {

    row <- c(na.omit(unlist(unname(NGT[i, ]))))

    if (length(row) >= (nfix + 1)) {
      rslt <- rbind(rslt, row[1 : (nfix + 1)])
    } else {
      break
    }
  }

  return(rslt)
}

# ------------------------------------------------------------------------------
# plot different stacking versions

# 1. stack_all vs stack_old vs stack_new

filter.window <- 11

ngt <- processAccumulationNGT() %>%
  dplyr::filter(Year >= 1502) %>%
  filterData(window = filter.window)

# stack_all
acc1 <- ngt %>%
  stackAllCores()
# stack_new
acc2 <- data.frame(
  Year = ngt$Year,
  stack = rowMeans(ngt[, 2 : 6], na.rm = TRUE))
# stack_new w/o B21-12
acc3 <- data.frame(
  Year = ngt$Year,
  stack = rowMeans(ngt[, c(2, 4 : 6)], na.rm = TRUE))
# stack_old
acc4 <- data.frame(
  Year = ngt$Year,
  stack = rowMeans(ngt[, 7 : 12], na.rm = TRUE))
# frozen stack (n=4, from 2010, w/o B21-12)
tmp <- ngt %>%
  dplyr::filter(Year <= 2010) %>%
  dplyr::select(-`B21_12`) %>%
  selectFixedNumber(nfix = 4)
acc5 <- data.frame(Year = tmp[, 1], stack = rowMeans(tmp[, -1]))

cols <- c("black", "#1b9e77", "#d95f02", "#7570b3")

grfxtools::Quartz(height = 6, width = 10,
                  file = "./zzz/ngt_acc_stack_versions_zoom.pdf")

plot(acc1, type = "n",
     xlim = c(1850, 2020), ylim = c(-10, 60), xlab = "", ylab = "")

lines(acc1, lwd = 3, col = cols[1])
lines(acc2, lwd = 1.5, col = cols[2])
lines(acc3, lwd = 2.5, col = "grey", lty = 2)
lines(acc4, lwd = 1.5, col = cols[3])
lines(acc5, lwd = 1.5, col = cols[4])

legend("topleft",
       c("stack_all", "stack_new_only", "stack_old_only",
         "frozen stack (n = 4, until 2010 CE, w/o B21-12)",
         "stack_new_only w/o B21-12"),
       lwd = 2, col = c(cols, "grey"), lty = c(rep(1, 4), 2), bty = "n")

dev.off()

# ------------------------------------------------------------------------------
# plot singe records

filter.window <- 11

ngt <- processAccumulationNGT() %>%
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

## legend("bottomleft", col = 1 : 11, lwd = 1.5, lty = 1, bty = "n",
##        legend = c("B18-2012", "B21-2012", "B23-2012", "B26-2012", "NGRIP-2012",
##                   "NEEM", "B16", "B18", "B21", "B26", "B29"))

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
