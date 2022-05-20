##
## aim:
## compare the NGT stack with Greenlandic Blocking Index data.
## relation:
## NGT paper supplementary; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

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

annualGBI <- readGBI() %>%
  filterData(window = filter.window) %>%
  dplyr::select("Year", "annual")

# ------------------------------------------------------------------------------
# scatter plots

permil2temperature <- 1 / 0.67

x1 <- subsetData(annualGBI, 2011 : 1851, "annual")
y1 <- subsetData(filteredStackedNGT, 2011 : 1851, "stack") * permil2temperature

x2 <- subsetData(annualGBI, 2011 : 1871, "annual")
y2 <- subsetData(filteredMAR, 2011 : 1871, "melt")

x3 <- subsetData(filteredStackedNGT, 2011 : 1871, "stack") * permil2temperature
y3 <- subsetData(filteredMAR, 2011 : 1871, "melt")

lab1 <- "Annual mean GBI"
lab2 <- grfxtools::LabelAxis("NGT-2012 temperature anomaly", unit = "celsius")
lab3 <- grfxtools::LabelAxis("Melt runoff anomaly", unit = "Gt",
                             unit.type = "trend", time.unit = "yr")

grfxtools::Quartz(file = "fig/supplement-gbi-comparison.pdf",
                  height = 4, width = 12)
par(mfrow = c(1, 3))

plot(x1, y1, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = c(-0.5, 1), ylim = c(-1.5, 3))

axis(1)
axis(2)
mtext(lab1, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(lab2, side = 2, line = 3.5, cex = par()$cex.lab * par()$cex, las = 0)
mtext("a", side = 3, adj = 0.02, padj = 0.5,
      line = -1.5, font = 2, cex = par()$cex.lab * par()$cex)

points(x1, y1, pch = 19)

plot(x2, y2, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = c(-0.5, 1), ylim = c(-100, 150))

axis(1)
axis(2)
mtext(lab1, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(lab3, side = 2, line = 3.5, cex = par()$cex.lab * par()$cex, las = 0)
mtext("b", side = 3, adj = 0.02, padj = 0.5,
      line = -1.5, font = 2, cex = par()$cex.lab * par()$cex)

points(x2, y2, pch = 19)

plot(x3, y3, type = "n", axes = FALSE,
     xlab = "", ylab = "", xlim = c(-1.5, 3), ylim = c(-100, 150))

axis(1)
axis(2)
mtext(lab2, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
mtext(lab3, side = 2, line = 3.5, cex = par()$cex.lab * par()$cex, las = 0)
mtext("c", side = 3, adj = 0.02, padj = 0.5,
      line = -1.5, font = 2, cex = par()$cex.lab * par()$cex)

points(x3, y3, pch = 19)
lines(x <- seq(-1, 3, 0.1), coef(lm(y3 ~ x3 + 0)) * x,
      col = "dodgerblue4", lwd = 2)

dev.off()

