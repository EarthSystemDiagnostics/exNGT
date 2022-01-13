##
## aim:
## plot point correlations between HadCrut/20CR fields and NGT-2012/Arctic2k
## relation:
## NGT paper supplementary; https://github.com/EarthSystemDiagnostics/exNGT
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

    theme(legend.key.height = unit(0.75, units = "inches"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          text = element_text(size = 18))

  p <- grfxtools::ggpolar(pole = "N", max.lat = 90, min.lat = 50,
                          n.lat.labels = 4,
                          longitude.spacing = 45,
                          land.fill.colour = "transparent",
                          size.outer = 0.5,
                          lat.ax.labs.pos = 180, ax.labs.size = 4.5,
                          data.layer = p)

  return(p)

}

# extract summar values for our ice core region on Greenland
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
# Plot 20CR results

# data
dat <- readRDS("out/supplement-point-correlations-20CR.rds")
dat <- dat[c(3, 4, 1, 2)]

# summary statistics for Greenland region
extractGreenlandRegion(dat$ngt.ann)
extractGreenlandRegion(dat$ngt.11y)

# plot panel labels
labels <- c(expression("(" * bold("a") * ") " * "NGT-2012, 11-yr r.m."),
            expression("(" * bold("b") * ") " * "Arctic2k, 11-yr r.m."),
            expression("(" * bold("c") * ") " * "NGT-2012, annual"),
            expression("(" * bold("d") * ") " * "Arctic2k, annual"))

# create plots
ggplt <- lapply(dat, plotMap, markNonsignificance = TRUE)

grfxtools::Quartz(height = 12, width = 14)
egg::ggarrange(plots = ggplt, nrow = 2, ncol = 2, labels = labels,
               label.args = list(gp = grid::gpar(cex = 1.25)))

filename <- "supplement_point_cor_20cr.pdf"
dev.copy2pdf(file = file.path("fig", filename))

