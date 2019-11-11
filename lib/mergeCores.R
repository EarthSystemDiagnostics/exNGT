doMerge <- function(site = "B18", data, method = 1, adjustMean = FALSE) {

  if (!method %in% c(1, 2))
    stop("Only 'method = 1' or 'method = 2' available.")

  startYear <- data$Year[1]
  cores <- calculateOverlapStatistics(data = data, site = site)

  if (adjustMean) {  
    cores$dat$pairData[, 3] <-
      cores$dat$pairData[, 3] - cores$stat$diffMeanOverlap
  }

  if (method == 1) {

    i <- match(startYear : cores$stat$startOverlap, data$Year)
    merged <- c(cores$dat$pairData[i, 3], cores$dat$pairData[-i, 2])

  } else {

    i <- match(startYear : (cores$stat$endOverlap + 1), data$Year)
    merged <- c(cores$dat$pairData[i, 3], cores$dat$pairData[-i, 2])

  }

  return(merged)

}

mergeCores <- function(data, sites = c("B18", "B21", "B23", "B26"),
                       method = 1, adjustMean = FALSE) {

  if (!method %in% c(1, 2))
    stop("Only 'method = 1' or 'method = 2' available.")

  merged <- sapply(sites, doMerge, data, method, adjustMean)
  merged <- as.data.frame(cbind(Year = data$Year, merged))

  return(merged)

}

