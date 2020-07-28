#' Stack cores
#'
#' Average (stack) the oxygen isotope data of the specified firn cores.
#'
#' @param data data frame with the relevant oxygen isotope records.
#' @param sites vector of character IDs of the sites for which the oxygen
#'   isotope records shall be averaged. The specified sites must be included in
#'   the given data.
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'   stripped before the computation proceeds; defaults to \code{TRUE}.
#' @return a data frame with two named columns 'Year' and 'stack' with the age
#'   scale and the averaged data.
#' @author Thomas Münch
#'
stackCores <- function(data, sites, na.rm = TRUE) {

  if (any(is.na(match(sites, names(data)))))
    stop("Requested site(s) not found in given data.")

  data.frame(Year = data$Year,
             stack = rowMeans(data[, sites], na.rm = na.rm))

}

#' Stack old and new cores
#'
#' Average the oxygen isotope records of (1) the redrilled cores and (2) all
#' available old cores.
#'
#' Note that the isotope records from the NEGIS and NEEM site are special in the
#' sense that they are both neither "true" nor "old" records; therefore they can
#' be excluded from the stack of all records.
#'
#' @param data data frame with the NGT and associated oxygen isotope records.
#' @param use_NEGIS_NEEM logical; indicate whether to include the records from
#'   the NEGIS and NEEM sites in the stack; defaults to \code{TRUE}.
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'   stripped before the computation proceeds; defaults to \code{TRUE}.
#' @return a data frame with three named columns 'Year', 'stack' and 'stack_12'
#'   with the age scale and the data from averaging the old and new cores,
#'   respectively.
#' @author Thomas Münch
#'
stackOldAndNew <- function(data, use_NEGIS_NEEM = TRUE, na.rm = TRUE) {

  sitesNew <- c("B18_12", "B21_12", "B23_12", "B26_12", "NGRIP_12")

  if (use_NEGIS_NEEM) {
    sitesOld <- names(data)[-match(c("Year", sitesNew),
                                   names(data))]
  } else {
    sitesOld <- names(data)[-match(c("Year", sitesNew, "NEGIS", "NEEM"),
                                   names(data))]
  }

  data.frame(
    Year     = data$Year,
    stack    = stackCores(data, sitesOld, na.rm = na.rm)$stack,
    stack_12 = stackCores(data, sitesNew, na.rm = na.rm)$stack
  )

}

#' Stack all cores
#'
#' Average all available oxygen isotope records of the NGT data set.
#'
#' Note that the isotope records from the NEGIS and NEEM site are special in the
#' sense that they are both neither "true" nor "old" records; therefore they can
#' be excluded from the stack of all records.
#'
#' @param data data frame with the NGT and associated oxygen isotope records.
#' @param use_NEGIS_NEEM logical; indicate whether to include the records from
#'   the NEGIS and NEEM sites in the stack; defaults to \code{TRUE}.
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'   stripped before the computation proceeds; defaults to \code{TRUE}.
#' @return a data frame with two named columns 'Year' and 'stack' with the age
#'   scale and the data from averaging all available cores.
#' @author Thomas Münch
#'
stackAllCores <- function(data, use_NEGIS_NEEM = TRUE, na.rm = TRUE) {

  sites <- names(data)[-1]

  if (!use_NEGIS_NEEM) {
    sites <- sites[-match(c("NEGIS", "NEEM"), sites)]
  }

  stackCores(data, sites, na.rm = na.rm)

}

#' Stack extended cores
#'
#' Average the set of oxygen isotope records created from merging the pairs of
#' old and new (re-drilled) cores from the same sites together with the other
#' available NGT and associated records.
#'
#' Note that the isotope records from the NEGIS and NEEM site are special in the
#' sense that they are both neither "true" nor "old" records; therefore they can
#' be excluded from the stack of all records.
#'
#' @param dataMerged data frame with a year column and the data columns from
#'   merging the pairs of old and new (re-drilled) cores from the same sites.
#' @param data original data frame with the NGT and associated oxygen isotope
#'   records.
#' @param use_NEGIS_NEEM logical; indicate whether to include the records from
#'   the NEGIS and NEEM sites in the stack; defaults to \code{TRUE}.
#' @param stack logical; whether to actually stack (average) the records.
#'   Defaults to \code{TRUE}. For \code{FALSE} the data frame with the records
#'   contributing to the stack is returned.
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'   stripped before the computation proceeds; defaults to \code{TRUE}.
#' @return a data frame with two named columns 'Year' and 'stack' with the age
#'   scale and the data from averaging the cores (for \code{stack = TRUE}).
#' @author Thomas Münch
#'
stackExtendedCores <- function(dataMerged, data, use_NEGIS_NEEM = TRUE,
                               stack = TRUE, na.rm = TRUE) {

  sites <- c("B18", "B21", "B23", "B26", "NGRIP")
  sites <- c(sites, paste0(sites, "_12"))

  if (!use_NEGIS_NEEM) {sites <- c(sites, "NEGIS", "NEEM")}

  data[, sites] <- NULL
  dataMerged$Year <- NULL

  data <- cbind(data, dataMerged)

  if (stack) {
    stackCores(data, names(data)[-1], na.rm = na.rm)
  } else {
    data
  }

}

#' Produce main NGT stack
#'
#' Wrapper function to produce the NGT stack used for the main part of the
#' paper.
#'
#' @param ngt data frame with the original or filtered NGT data.
#' @param stack logical; whether to actually stack (average) the records.
#'   Defaults to \code{TRUE}. For \code{FALSE} the data frame with the records
#'   contributing to the stack is returned.
#' @return data frame with the NGT age scale as the first and the stacked
#'   isotope values as the second column (for \code{stack = TRUE}).
#' @author Thomas Münch
#'
stackNGT <- function(ngt, stack = TRUE) {

  ngt %>%
    mergeCores(adjustMean = TRUE, method = 1) %>%
    stackExtendedCores(ngt, stack = stack)

}
