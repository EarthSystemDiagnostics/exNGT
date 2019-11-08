#Function to derive statistics in overlap period of extended NGT records
#by Maria Hoerhold and Thomas Muench
#' @title  deliver statistics of overlap
#' @description read and compare old and extended records of site 
#' @param Input NGT dataframe
#' @param Bname/ Site = name of site
#' @return statistics of overlap period
#' 

# function to calculate the standard error
sd_error<-function(x)
{result <- sd(x)/ sqrt(length(x))
return(result)}

calculateOverlapStatistics<-function(data, site = "B18") {
  pairID <- c("Year", site, paste(site, "12", sep = "_"))
  pairData <- data[, pairID]
  overlap <- na.omit(pairData)
  
  output <- list(
    startOverlap = min(overlap$Year),
    endOverlap = max(overlap$Year),
    meanOverlap = colMeans(overlap[, -1]),
    diffMeanOverlap = diff(colMeans(overlap[, -1])),
    corrMeanOverlap = cor(overlap[, -1]),
    sderrorOverlap = apply(overlap[, -1], 2, sd_error))
  
  return(output)
  
}




