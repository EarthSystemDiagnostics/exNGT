##
## aim:
## initialize the "exNGT" R function bundle by loading all needed
## packages and library functions.
## relation:
## https://github.com/EarthSystemDiagnostics/exNGT
##
## Maria Hoerhold and Thomas Muench, 11/2019
##


# ------------------------------------------------------------------------------
# Load all required packages

## none to load so far


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
# Source the "exNGT" library directory

sourceDir("lib")

