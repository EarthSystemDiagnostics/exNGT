##
## aim:
## analyse overlap period of existing and re-drilled cores.
## relation:
## https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

library(dplyr)

# ------------------------------------------------------------------------------
# Get relevant data

# all NGT anomaly records
NGT <- processNGT()

# add stacks of old and new records
NGT <- cbind(NGT, stackOldAndNew(NGT)[-1])

# ------------------------------------------------------------------------------
# Loop over running mean window sizes and calculate overlap stats for each core

# set running mean window sizes
w <- c(1, 3, 5, 7, 11, 21)

# record IDs which have old and new records available
sites <- c("B18", "B21", "B23", "B26", "NGRIP", "stack")

# calculate running mean averages
filteredRecords <- sapply(w, simplify = FALSE, function(x) {
  filterData(NGT, window = x)
})

# calculate overlap statistics for each core and window size
overlapStats <- sapply(sites, simplify = FALSE, FUN = function(s) {

  res <- sapply(w, simplify = FALSE, function(x) {

    rec <- filteredRecords[[match(x, w)]]

    res <- data.frame(window = x)
    res <- cbind(res, calculateOverlapStatistics(rec, site = s)$stat)

    return(res)

  })

  do.call(rbind, res)
})

# arrange as single data frame
overlapStats <- do.call(function(...){rbind(..., make.row.names = FALSE)},
                        overlapStats)

# ------------------------------------------------------------------------------
# Show certain variables

# show all overlap correlations
overlapStats %>%
  select(corePair, window, corrOverlap)

# show overlap correlations only for certain window size
overlapStats %>%
  filter(window == 11) %>%
  select(corePair, corrOverlap)
  
