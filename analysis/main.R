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

asp <- 7 / 15
w <- 18.3
natfig(file = "./fig/main-figure01ac.pdf",
       height = asp * w, width = w)

makeFigure01()
dev.off()

asp <- 1
w <- 0.3 * 18.3
natfig(file = "./fig/main-figure01b-raw.pdf",
       height = asp * w, width = w)

plotMap()
dev.off()

# ------------------------------------------------------------------------------
# Figure 03 - spectrum and coherence

asp <- 7 / 8
w <- 8.9
natfig(file = "./fig/main-figure03.pdf",
       height = asp * w, width = w, mar = c(5, 5, 0.5, 5))

makeFigure03()
dev.off()

# ------------------------------------------------------------------------------
# Figure 04a - histogram

asp <- 6.2 / 7.5
w <- 8.6
natfig(file = "./fig/main-figure04a.pdf",
       height = asp * w, width = w)

makeFigure04()
dev.off()

# ------------------------------------------------------------------------------
# Run additional analyses

rmarkdown::render(input = "./analysis/analysis-overview.Rmd",
                  output_dir = "./out")
rmarkdown::render(input = "./analysis/probability-analysis.Rmd",
                  output_dir = "./out")
rmarkdown::render(input = "./analysis/overlap-analysis.Rmd",
                  output_dir = "./out")
