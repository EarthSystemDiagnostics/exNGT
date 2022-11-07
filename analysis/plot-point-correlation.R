##
## aim:
## plot point correlation between 20CR reanalysis field and NGT-2012/Arctic2k
## relation:
## NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

library(ggplot2)

# ------------------------------------------------------------------------------
# Helper functions

# polar stereographic plot
plotMap <- function(map, MAX = 1, MIN = -MAX, markNonsignificance = FALSE,
                    pval = 0.05, na.fill = "grey") {

  colour.scale <- grfxtools::ColorPal("RdBu", rev = TRUE)

  cor.min <- MIN - 0.01
  cor.max <- MAX + 0.01

  map$r[map$r > MAX] <- MAX
  map$r[map$r < MIN] <- MIN

  if (markNonsignificance) {
    map$r[map$p > pval] <- NA
  }

  p <- ggplot() +

    geom_tile(aes(x = lon, y = lat, fill = r),
              data = map, colour = "transparent") +

    scale_fill_gradientn(colours = colour.scale,
                         na.value = na.fill,
                         limits = c(cor.min, cor.max),
                         name = "Correlation") +

    theme(legend.key.height = unit(0.4, units = "inches"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))

  p <- grfxtools::ggpolar(pole = "N", max.lat = 90, min.lat = 49.5,
                          lat.ax.vals = c(60, 70, 80),
                          longitude.spacing = 45,
                          land.fill.colour = "transparent",
                          size = pt2mm(1.338 * 0.495), ax.labs.size = pt2mm(6),
                          size.outer = pt2mm(1.338 * 0.495),
                          size.axes = pt2mm(1.32 * 0.33),
                          lat.ax.labs.pos = 180, data.layer = p)

  return(p)

}

# extract summary values for our ice core region on Greenland
extractGreenlandRegion <- function(x) {

  corePos <- loadPositions() %>%
    dplyr::filter(Identifier == "core")

  minLat <- x$lat[which.min(abs(min(corePos$Latitude) - x$lat))]
  maxLat <- x$lat[which.min(abs(max(corePos$Latitude) - x$lat))]
  minLon <- x$lon[which.min(abs(min(corePos$Longitude) - x$lon))]
  maxLon <- x$lon[which.min(abs(max(corePos$Longitude) - x$lon))]

  x %>%
    dplyr::filter(lat >= minLat & lat <= maxLat) %>%
    dplyr::filter(lon >= minLon & lon <= maxLon) %>%
    dplyr::summarise(rMean = mean(r), rMin = min(r), rMax = max(r), p = max(p))

}

# ------------------------------------------------------------------------------
# Plot 20CR main results with filtered data

# data
dat <- readRDS("out/point-correlations-20CR.rds")
dat <- dat[c(3, 4, 1, 2)]

# summary statistics for Greenland region
extractGreenlandRegion(dat$ngt.ann)
extractGreenlandRegion(dat$ngt.11y)

# correlation with area-weighted mean reanalysis time series
filter.window <- 11
stackedNGT <- processNGT() %>%
  stackNGT()
filteredStackedNGT <- processNGT() %>%
  filterData(window = filter.window) %>%
  stackNGT()
tcr <- readTwenCR() %>%
  filterData(window = filter.window)

res <- estimateCorrelation(stackedNGT, filteredStackedNGT, tcr,
                           filter.window = filter.window,
                           analysis.period = 2011 : 1836, nmc = 10000)
sprintf("r = %1.2f (p = %1.4f)", res$r, res$p)

# plot panel labels
labels <- c(expression(bold("a")), expression(bold("b")))

# create plots
ggplt1 <- lapply(dat[1 : 2], plotMap, markNonsignificance = TRUE)
ggplt2 <- lapply(dat[3 : 4], plotMap, markNonsignificance = TRUE)

asp <- 6 / 14
w <- 18.3
natfig(height = asp * w, width = w)

egg::ggarrange(plots = ggplt1, nrow = 1, ncol = 2, labels = labels,
               label.args = list(gp = grid::gpar(fontsize = 8)))

filename <- "main-figure02.pdf"
dev.copy2pdf(file = file.path("fig", filename))

# ------------------------------------------------------------------------------
# Plot 20CR supplementary results with annual data

dat <- readRDS("out/point-correlations-20CR-with-20CR.rds")

labels <- c(expression(bold("a")), expression(bold("b")),
            expression(bold("c")), expression(bold("d")))

# create plots
ggplt3 <- lapply(dat, plotMap, markNonsignificance = TRUE)

asp <- 12 / 14
w <- 18.3
natfig(height = asp * w, width = w)

egg::ggarrange(plots = c(ggplt2, ggplt3), nrow = 2, ncol = 2, labels = labels,
               label.args = list(gp = grid::gpar(fontsize = 8)))

filename <- "supplement_point_cor_20cr_annual.pdf"
dev.copy2pdf(file = file.path("fig", filename))
