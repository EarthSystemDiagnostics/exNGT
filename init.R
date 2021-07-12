##
## aim:
## initialize the "exNGT" R function bundle by loading all needed
## packages and library functions.
## relation:
## https://github.com/EarthSystemDiagnostics/exNGT
##
## Maria Hoerhold and Thomas Muench, 05/2021
##


# ------------------------------------------------------------------------------
# Check if required packages are available

required.packages <- c("grfxtools", "proxysnr", "tibble", "ggplot2", "magrittr",
                       "dplyr", "rmarkdown", "zoo")

for (pkg in required.packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package \"%s\" needed for \"exNGT\". ", pkg),
         "Please install it using \"dependencies.R\".", call. = FALSE)
  }
}

# ------------------------------------------------------------------------------
# Define function to source entire directory with .R files

sourceDir <- function (path, trace = TRUE, local = FALSE, ...) {
  cat("Sourcing files...\n")
  for (nm in list.files(path, pattern = "[.][Rr]$")) {
    if (trace) cat(nm, ":")
    source(file.path(path, nm), local = local, ...)
    if (trace) cat("\n")
  }
}

# ------------------------------------------------------------------------------
# Attach necessary packages and source the "exNGT" library directory

library("magrittr")
sourceDir("lib")

