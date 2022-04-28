readGBI <- function(path = "data/gbi.monthly.csv") {

  read.csv(file = path, skip = 4) %>%
    dplyr::arrange(dplyr::desc(Year))

}

filter.window <- 11

stackedNGT <- processNGT() %>%
  stackNGT()

filteredStackedNGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()

MAR <- readMAR() %>%
  dplyr::select("Year", "melt")

filteredMAR <- MAR %>%
  filterData(window = filter.window)

Arctic2k <- readArctic2k() %>%
  extendWithHadCrut() %>%
  dplyr::select("Year", "TempAnomaly")

filteredArctic2k <- Arctic2k %>%
  filterData(window = filter.window)

GBI <- readGBI() %>%
  filterData(window = filter.window)

annualGBI <- dplyr::select(GBI, "Year", "annual")
summerGBI <- dplyr::select(GBI, "Year", "summer")


estimateCorrelation(stackedNGT, filteredStackedNGT, annualGBI,
                    filter.window = filter.window,
                    analysis.period = 2011 : 1961)
estimateCorrelation(stackedNGT, filteredStackedNGT, annualGBI,
                    filter.window = filter.window,
                    analysis.period = 1960 : 1851)

estimateCorrelation(Arctic2k, filteredArctic2k, annualGBI,
                    filter.window = filter.window,
                    analysis.period = 2011 : 1961)
estimateCorrelation(Arctic2k, filteredArctic2k, annualGBI,
                    filter.window = filter.window,
                    analysis.period = 1960 : 1851)

estimateCorrelation(MAR, filteredMAR, summerGBI,
                    filter.window = filter.window,
                    analysis.period = 2011 : 1961)
estimateCorrelation(MAR, filteredMAR, summerGBI,
                    filter.window = filter.window,
                    analysis.period = 1960 : 1871)

# ------------------------------------------------------------------------------
# scatter plots

t1 <- 2011 : 1961
t2 <- 1960 : 1851
t3 <- 1960 : 1871

x1 <- subsetData(filteredStackedNGT, t1, "stack")
y1 <- subsetData(annualGBI, t1, "annual")
x2 <- subsetData(filteredStackedNGT, t2, "stack")
y2 <- subsetData(annualGBI, t2, "annual")

x3 <- subsetData(filteredArctic2k, t1, "TempAnomaly")
y3 <- subsetData(annualGBI, t1, "annual")
x4 <- subsetData(filteredArctic2k, t2, "TempAnomaly")
y4 <- subsetData(annualGBI, t2, "annual")

x5 <- subsetData(filteredMAR, t1, "melt")
y5 <- subsetData(summerGBI, t1, "summer")
x6 <- subsetData(filteredMAR, t3, "melt")
y6 <- subsetData(summerGBI, t3, "summer")

col1 <- "firebrick3"
col2 <- "black"

xlab1 <- grfxtools::LabelAxis("NGT-2012")
xlab2 <- grfxtools::LabelAxis("Arctic 2k", unit = "celsius")
xlab3 <- grfxtools::LabelAxis("Melt runoff", unit = "Gt",
                              unit.type = "trend", time.unit = "yr")
ylab1 <- "Annual mean GBI"
ylab2 <- "Summer GBI"

grfxtools::Quartz(file = "fig/supplement-gbi-correlation.pdf",
                  height = 6, width = 18)
par(mfrow = c(1, 3))

plot(x1, y1, type = "n", xlab = xlab1, ylab = ylab1,
     xlim = c(-1, 2), ylim = c(-0.5, 1))

points(x1, y1, pch = 19, col = col1)
points(x2, y2, pch = 19, col = col2)

legend("bottomright", c("1961 - 2011 CE", "1851 - 1960 CE"),
       col = c(col1, col2), pch = 19, bty = "n")

plot(x3, y3, type = "n", xlab = xlab2, ylab = ylab1,
     xlim = c(-1.5, 2), ylim = c(-0.5, 1))

points(x3, y3, pch = 19, col = col1)
points(x4, y4, pch = 19, col = col2)

legend("bottomright", c("1961 - 2011 CE", "1851 - 1960 CE"),
       col = c(col1, col2), pch = 19, bty = "n")

plot(x5, y5, type = "n", xlab = xlab3, ylab = ylab2,
     xlim = c(-75, 150), ylim = c(-0.5, 1))

points(x5, y5, pch = 19, col = col1)
points(x6, y6, pch = 19, col = col2)

legend("bottomright", c("1961 - 2011 CE", "1871 - 1960 CE"),
       col = c(col1, col2), pch = 19, bty = "n")

dev.off()
