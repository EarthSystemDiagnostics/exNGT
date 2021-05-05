##
## aim:
## install all package dependencies for "exNGT"
## relation:
## https://github.com/EarthSystemDiagnostics/exNGT
##
## Maria Hoerhold and Thomas Muench, 05/2021
##

# ------------------------------------------------------------------------------
# install dependencies for all analyses and plotting routines

install.packages("dplyr")
install.packages("ggplot2")
install.packages("tibble")

# install.packages("remotes")
remotes::install_github("EarthSystemDiagnostics/proxysnr")
remotes::install_github("EarthSystemDiagnostics/grfxtools")

# ------------------------------------------------------------------------------
# install optional dependencies needed only for processing of the raw data

install.packages("pangaear")
install.packages("pdftools")
install.packages("plyr")
install.packages("purrr")
install.packages("readr")
install.packages("stringr")
