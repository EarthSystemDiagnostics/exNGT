#' Merge old and new records
#'
#' Merge the old and new (re-drilled) oxygen isotope records of a specific site,
#' either at the start or the end of the overlap period and with or without
#' adjusting the mean values of the records before merging.
#'
#' @param site character string of the site ID for which to merge the old and
#'   new records; must match the list of sites for which re-drilled records are
#'   available.
#' @param data data frame with the NGT and associated oxygen isotope records.
#' @param mergePoint character signalling the time point of merging: for
#'   \code{"start"} (the default) the new record is merged with the old record
#'   at the start (i.e. oldest age) of their common overlap period, for
#'   \code{"end"} the records are merged at the end (i.e. one year after the
#'   youngest age) of the overlap period.
#' @param adjustMean logical; shall the mean values of the two records within
#'   the common overlap period be adjusted prior to the merging? Defaults to
#'   \code{FALSE}. If \code{TRUE}, mean adjustment is always performed for the
#'   new record relative to the old one.
#' @return numeric vector with the merged record.
#' @author Thomas Münch
doMerge <- function(site = c("B18", "B21", "B23", "B26", "NGRIP", "stack"),
                    data, mergePoint = c("start", "end"), adjustMean = FALSE) {

  site <- match.arg(site)
  mergePoint <- match.arg(mergePoint)

  cores <- calculateOverlapStatistics(data = data, site = site)

  startYear <- data$Year[1]
  endYear   <- ifelse(mergePoint == "start", cores$stat$startOverlap,
                      cores$stat$endOverlap + 1)

  if (endYear > startYear) {
    stop("Cannot define merging point: ",
         "old and new records have same final age.", call. = FALSE)
  }
  if (endYear == dplyr::last(data$Year)) {
   stop("Cannot define merging point: ",
         "old and new records have same start age.", call. = FALSE)
  }

  if (adjustMean) {  
    cores$dat$pairData[, 3] <-
      cores$dat$pairData[, 3] - cores$stat$diffMeanOverlap
  }

  i <- match(startYear : endYear, data$Year)
  merged <- c(cores$dat$pairData[i, 3], cores$dat$pairData[-i, 2])

  return(merged)

}

#' Merge a set of old and new records
#'
#' Wrapper function to merge all suitable pairs of old and new oxygen isotope
#' records from the NGT data set.
#'
#' @param data data frame with the NGT and associated oxygen isotope records.
#' @param sites vector of the character IDs of the sites for which the merging
#'   shall be conducted; defaults to the list of sites with available paired
#'   records.
#' @param mergePoint character signalling the time point of merging: for
#'   \code{"start"} (the default) the new record is merged with the old record
#'   at the start (i.e. oldest age) of their common overlap period, for
#'   \code{"end"} the records are merged at the end (i.e. one year after the
#'   youngest age) of the overlap period.
#' @param adjustMean logical; shall the mean values of the two records within
#'   the common overlap period be adjusted prior to the merging? Defaults to
#'   \code{FALSE}. If \code{TRUE}, mean adjustment is always performed for the
#'   new record relative to the old one.
#' @return a data frame with an age column and the merged versions of the
#'   records specified in \code{sites}.
#' @author Thomas Münch
mergeCores <- function(data,
                       sites = c("B18", "B21", "B23", "B26", "NGRIP"),
                       mergePoint = c("start", "end"), adjustMean = FALSE) {

  mergePoint <- match.arg(mergePoint)

  if (any(is.na(match(sites, names(data)))))
    stop("Requested site(s) not found in given data.")

  merged <- sapply(sites, doMerge, data, mergePoint, adjustMean)
  merged <- as.data.frame(cbind(Year = data$Year, merged))

  return(merged)

}
