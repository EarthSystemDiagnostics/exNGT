#' Plot histogram
#'
#' Plot a histogram of isotope values marking certain time periods in the
#' isotope record.
#'
#' @param x a data frame with two columns with the time scale as the first and
#'   the isotope values as the second column.
#' @param analysis.period optional vector of time points to subset the data over
#'   which the full histogram shall be calculated. Set to \code{NULL} to use the
#'   full range of \code{x}.
#' @param p1 vector of time points for marking a first record period.
#' @param p2 vector of time points for marking a second record period.
#' @param p3 vector of time points for marking a third record period.
#' @param p4 vector of time points for marking a fourth record period.
#' @param breaks a vector of break points between the histogram cells.
#' @param xlim x axis range; defaults to the range of \code{breaks}.
#' @param ylim y axis range.
#' @param col vector of four distinct colours to mark the four time periods in
#'   the histogram.
#' @param xlab the x axis label; setting it to "anomaly" or "slope" produces
#'   isotope-specific labels for isotope anomaly or isotope slope histograms.
#' @param ylab the y axis label.
#' @param alpha opacity value within [0,1] for the histogram \code{colours}.
#' @param plot.quantiles logical; set to \code{TRUE} to plot vertical lines for
#'   the 5, 50 and 95 % quantiles of the full histogram.
#' @author Thomas Münch
#'
plotHistogram <- function(x, analysis.period = 1000 : 2011,
                          p1 = 2011 : 1997, p2 = 1938 : 1924,
                          p3 = 1884 : 1870, p4 = 1424 : 1410,
                          breaks = seq(-2.5, 2.5, 0.25),
                          xlim = range(breaks), ylim = c(0, 0.3),
                          col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"),
                          xlab = "anomaly", ylab = "Relative counts",
                          alpha = 0.6, plot.quantiles = TRUE) {

  if (length(col) != 4) stop("Supply four distinct colours.")

  if (xlab == "anomaly") {
    xlab = bquote(delta^{"18"} * "O anomaly (\u2030)")
  } else if (xlab == "slope") {
    xlab = bquote(delta^{"18"} * "O slope (\u2030 " * "yr"^{-1} * ")")
  }

  if (!is.null(analysis.period)) {
    x <- x[match(analysis.period, x[, 1]), ]
  }

  col <- adjustcolor(col, alpha = alpha)

  nobs <- length(na.omit(x[, 2]))

  x0 <- x[, 2]

  x1 <- x[match(p1, x[, 1]), 2]
  x2 <- x[match(p2, x[, 1]), 2]
  x3 <- x[match(p3, x[, 1]), 2]
  x4 <- x[match(p4, x[, 1]), 2]

  h0 <- hist(x0, breaks = breaks, plot = FALSE)
  h1 <- hist(x1, breaks = breaks, plot = FALSE)
  h2 <- hist(x2, breaks = breaks, plot = FALSE)
  h3 <- hist(x3, breaks = breaks, plot = FALSE)
  h4 <- hist(x4, breaks = breaks, plot = FALSE)

  h0$counts <- h0$counts / nobs
  h1$counts <- h1$counts / nobs
  h2$counts <- h2$counts / nobs
  h3$counts <- h3$counts / nobs
  h4$counts <- h4$counts / nobs

  plot(h0, xlim = xlim, ylim = ylim,
       main = "", xlab = "", ylab = "",
       col = adjustcolor("black", alpha = 0.2), yaxs = "i")

  mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
  mtext(ylab, side = 2, line = 3.5, cex = par()$cex.lab * par()$cex, las = 0)

  plot(h1, xlim = xlim, ylim = ylim,
       main = "", xlab = "", ylab = "",
       col = col[1], add = TRUE)

  plot(h2, xlim = xlim, ylim = ylim,
       main = "", xlab = "", ylab = "",
       col = col[2], add = TRUE)

  plot(h3, xlim = xlim, ylim = ylim,
       main = "", xlab = "", ylab = "",
       col = col[3], add = TRUE)

  plot(h4, xlim = xlim, ylim = ylim,
       main = "", xlab = "", ylab = "",
       col = col[4], add = TRUE)

  if (plot.quantiles) {
    abline(v = quantile(x0, probs = c(0.05, 0.5, 0.95), na.rm = TRUE),
           col = "black", lty = 5, lwd = 1)
  }

  lab0 <- "Full data"
  lab1 <- sprintf("%s to %s", min(p1), max(p1))
  lab2 <- sprintf("%s to %s", min(p2), max(p2))
  lab3 <- sprintf("%s to %s", min(p3), max(p3))
  lab4 <- sprintf("%s to %s", min(p4), max(p4))

  legend("topright", c(lab0, lab1, lab2, lab3, lab4),
         col = c(adjustcolor(1, 0.2), col), lty = 1, lwd = 10, bty = "n")

}
