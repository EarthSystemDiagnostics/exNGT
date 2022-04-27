##
## aim:
## make analyses and plots for the main paper part.
## relation:
## NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

# ------------------------------------------------------------------------------
# Figure 01 - time series, spectra and map

grfxtools::Quartz(file = "./fig/main-figure01.pdf",
                  height = 7.5, width = 15)

makeFigure01()
dev.off()

grfxtools::Quartz(file = "./fig/supplement-map-raw.pdf",
                  height = 6, width = 6)

plotMap()
dev.off()

# ------------------------------------------------------------------------------
# Figure 02 - histogram

grfxtools::Quartz(file = "./fig/main-figure02.pdf",
                  height = 6.2, width = 7.5)

makeFigure02()
dev.off()

# ------------------------------------------------------------------------------
# Run additional analyses

rmarkdown::render(input = "./analysis/analysis-overview.Rmd",
                  output_dir = "./out")
rmarkdown::render(input = "./analysis/probability-analysis.Rmd",
                  output_dir = "./out")
rmarkdown::render(input = "./analysis/overlap-analysis.Rmd",
                  output_dir = "./out")
