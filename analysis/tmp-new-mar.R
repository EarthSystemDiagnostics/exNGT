##
## aim:
## analyse new MAR data extendind till 1871.
## relation:
## NGT paper supplementary; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

# ------------------------------------------------------------------------------
# plot template

filter.window <- 11
permil2temperature.mid <- 1 / 0.67 # Greenland spatial slope
permil2temperature.low <- 1 / 1.1  # Masson-Delmotte et al. (2015)
permil2temperature.hig <- 2.1      # Vinther et al. (2009)

filteredStackedNGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()

filteredStackedTemperatureNGT.mid <- filteredStackedNGT %>%
  dplyr::mutate(stack = permil2temperature.mid * stack)
filteredStackedTemperatureNGT.low <- filteredStackedNGT %>%
  dplyr::mutate(stack = permil2temperature.low * stack)
filteredStackedTemperatureNGT.hig <- filteredStackedNGT %>%
  dplyr::mutate(stack = permil2temperature.hig * stack)

filteredMAR <- readMAR() %>%
  filterData(window = filter.window)  

col <- c("black", "#d95f02", "#7570b3")

grfxtools::Quartz()
par(mar = c(5, 5, 0.5, 5))

plot(filteredStackedTemperatureNGT.mid, type = "n", axes = FALSE,
     xlab = "", ylab = "", c(1871, 2011), ylim = c(-2.5, 3), xaxs = "i")

lines(filteredStackedTemperatureNGT.mid, lwd = 2, col = col[1])

axis(1)
axis(2)
mtext("Year CE", 1, 3.5, cex = par()$cex.lab)
mtext(grfxtools::LabelAxis("Temperature anomaly", unit = "celsius"),
      2, 3, cex = par()$cex.lab, las = 0)

grfxtools::Polyplot(filteredStackedTemperatureNGT.mid$Year,
                    filteredStackedTemperatureNGT.low$stack,
                    filteredStackedTemperatureNGT.hig$stack,
                    col = col[1], alpha = 0.2)
lines(filteredStackedTemperatureNGT.low, lwd = 1, col = col[1])
lines(filteredStackedTemperatureNGT.hig, lwd = 1, col = col[1])

lines(filteredMAR$Year, filteredMAR$t2m, lwd = 2, col = col[2])

par(new = TRUE)

plot(filteredMAR$Year, filteredMAR$melt, type = "l", axes = FALSE,
     col = col[3], lwd = 2,
     xlim = c(1871, 2011), ylim = c(-55, 300), xaxs = "i",
     xlab = "", ylab = "")

axis(4)
text(2030, 125, grfxtools::LabelAxis("Melt runoff anomaly", unit = "Gt",
                                     unit.type = "trend", time.unit = "yr"),
     srt = -90, xpd = NA, cex = par()$cex.lab * par()$cex)

legend("topleft",
       c("NGT-2012 temperature", "MAR3.5 temperature", "MAR3.5 melt runoff"),
       lty = 1, lwd = 2, col = col, bty = "n")

# ------------------------------------------------------------------------------
# NGT-MAR correlations

filter.window <- 11
permil2temperature <- 1 / 0.67

filteredMAR <- readMAR() %>%
  filterData(window = filter.window)

stackedTemperatureNGT <- processNGT() %>%
  stackNGT() %>%
  dplyr::mutate(stack = permil2temperature * stack)

filteredStackedTemperatureNGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT() %>%
  dplyr::mutate(stack = permil2temperature * stack)

analysis.period <- 2011 : 1871

estimateCorrelation(stackedTemperatureNGT,
                    filteredStackedTemperatureNGT,
                    filteredMAR[, c("Year", "t2m")],
                    filter.window = filter.window,
                    analysis.period = analysis.period,
                    nmc = 10000)
# <- r = 0.77 (p << 0.01)
estimateCorrelation(stackedTemperatureNGT,
                    filteredStackedTemperatureNGT,
                    filteredMAR[, c("Year", "melt")],
                    filter.window = filter.window,
                    analysis.period = analysis.period,
                    nmc = 10000)
# <- r = 0.63 (p = 0.004)

MAR <- readNewMAR() %>%
  makeAnomalies()
estimateCorrelation(MAR[, c("Year", "t2m")],
                    filteredMAR[, c("Year", "t2m")],
                    filteredMAR[, c("Year", "melt")],
                    filter.window = filter.window,
                    analysis.period = analysis.period,
                    nmc = 10000)
# <- r = 0.65 (p = 0.002)

