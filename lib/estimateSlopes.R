#' Estimate isotope slopes
#'
#' Estimate the linear trend slopes in an isotope record within a given
#' moving window.
#'
#' @param x a data frame with two columns with the time scale as the first and
#'   the isotope values as the second column.
#' @param window an integer value giving the size of the moving window within
#'   which the linear trend slopes are estimated; must be an odd number.
#' @return a data frame with two columns 'Year' and 'slope' with the time scale
#'   as the first and the slope values as the second column.
#' @author Thomas Münch
#'
estimateSlopes <- function(x, window = 21) {

  if (window < 3) stop("No sensible estimation window size.")
  if ((window %% 2) == 0) stop("Odd estimation window size expected.")
  
  slopes <- rep(NA, nrow(x))
  window.half <- (window - 1) / 2
  for (i in (window.half + 1) : (nrow(x) - window.half)) {

    j <- (i - window.half) : (i + window.half)

    slopes[i] <- coefficients(lm(x[j, 2] ~ x[j, 1]))[2]

  }

  res <- data.frame(Year = x[, 1], slope = slopes)
  return(res)

}
