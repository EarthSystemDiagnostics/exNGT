#Function to derive statistics in overlap period of extended NGT records
#by Maria Hoerhold and Thomas Muench
#' @title  deliver statistics of overlap
#' @description read and compare old and extended records of site 
#' @param Input NGT dataframe
#' @param Bname/ Site = name of site
#' @return statistics of overlap period
#' 

# function to calculate the standard error
sdError <- function(x, na.rm = FALSE) {
  sd(x, na.rm = na.rm) / sqrt(length(x))
}

calculateOverlapStatistics <- function(data, site = "B18") {

  if (!site %in% c("B18", "B21", "B23", "B26", "NGRIP", "stack"))
    stop("Unsuitable site requested.")

  pairID <- c("Year", site, paste(site, "12", sep = "_"))

  pairData <- data[, pairID]
  overlap <- na.omit(pairData)
  meanOverlap <- colMeans(overlap[, -1])
  
  output <- list(

    dat  = list(pairData = pairData, overlapData = overlap),
    stat = list(
      startOverlap = min(overlap$Year),
      endOverlap = max(overlap$Year),
      meanOverlap = meanOverlap,
      diffMeanOverlap = diff(meanOverlap),
      corrOverlap = cor(overlap[, 2], overlap[, 3]),
      sdErrorOverlap = apply(overlap[, -1], 2, sdError))
  )
  
  return(output)
  
}




