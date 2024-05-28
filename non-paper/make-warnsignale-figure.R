##
## aim:
##   create figures for chapter in Reklim book series "Warnsignale Klima":
##   - combine paper figures 1a (NGT-2012 part only) and 4a = book fig. 2;
##   - paper figure 1c as standalone version = book fig. 3a;
##   - all with German labels.
## relation:
##   NGT paper: <https://doi.org/10.1038/s41586-022-05517-z>
##   <https://github.com/EarthSystemDiagnostics/exNGT>
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

# ------------------------------------------------------------------------------
# Define figure functions

plotNGT <- function(filter.window = 11,
                    permil2temperature = 1 / 0.67) {

  if (filter.window < 1) {
    stop("Running mean filter window needs to be >= 1.", call. = FALSE)
  }

  if ((filter.window %% 2) == 0) {
    warning("Running mean filter window should be odd.", call. = FALSE)
  }

  # ----------------------------------------------------------------------------
  # Load time series data

  NGT <- processNGT()

  stackedNGT <- NGT %>%
    stackNGT() %>%
    dplyr::filter(Year >= 1000)

  filteredStackedNGT <- NGT %>%
    filterData(window = filter.window) %>%
    stackNGT() %>%
    dplyr::filter(Year >= 1000)

  # Linear regression data (annual means)

  t1 <- 1800 : 1000
  t2 <- 2011 : 1800

  regressionData <- list(
    data.frame(x = t1, y = subsetData(stackedNGT, t1, "stack")),
    data.frame(x = t2, y = subsetData(stackedNGT, t2, "stack"))
  )

  regressionModels <- regressionData %>%
    lapply(function(lst) {coef(lm(y ~ x, lst))})

  # ----------------------------------------------------------------------------
  # Make plots

  xlab  <- "Jahr n.d.Z."
  ylab.ngt.pm <- grfxtools::LabelAxis("NGT-2012")
  ylab.ngt.dc <- grfxtools::LabelAxis("NGT-2012", unit = "celsius")
  ylab.a2k <- grfxtools::LabelAxis("Arctic 2k", unit = "celsius")

  xlim <- c(1000, 2020)
  ylim.ngt <- c(-2.5, 2.5)

  xanml1 <- c(800, 2050)
  yanml <- rep(0, 2)

  x1 <- 845
  x2 <- 2180
  y1 <- 0.
  y2 <- 0.

  t.pi <- c(1000, 1800)
  t.rec <- c(2001, 2011)

  y.pi <- c(-1.85, 1.5)
  y.rec <- c(-1.85, 2.3)

  col <- c("black", "dodgerblue4", "wheat4")
  rect.col <- c("black", "orange2")

  plot(stackedNGT, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = xlim, ylim = ylim.ngt)

  rect(t.rec[1], y.rec[1], t.rec[2], y.rec[2],
       col = adjustcolor(rect.col[2], alpha = 0.8), border = NA)
  rect(t.pi[1], y.pi[1], t.pi[2], y.pi[2],
       col = adjustcolor(rect.col[1], alpha = 0.2), border = NA)

  lines(xanml1, yanml, lty = 2, lwd = Lwd(1.5), col = "darkgrey")

  lines(stackedNGT, col = "darkgrey")

  lines(filteredStackedNGT, col = col[1], lwd = Lwd(2.5))

  axisLwd(2, at = seq(-2, 2, 1))
  axisLwd(4, labels = seq(-3, 3, 1), at = seq(-3, 3, 1) / permil2temperature)

  text(x1, y1, ylab.ngt.pm, srt = +90, xpd = NA,
       cex = par()$cex.lab * par()$cex, col = col[1])
  text(x2, y2 / permil2temperature, ylab.ngt.dc, srt = -90, xpd = NA,
       cex = par()$cex.lab * par()$cex, col = col[1])

  lines(t1, regressionModels[[1]][1] + regressionModels[[1]][2] * t1,
        col = col[1], lwd = Lwd(2), lty = 2)
  lines(t2, regressionModels[[2]][1] + regressionModels[[2]][2] * t2,
        col = col[1], lwd = Lwd(2), lty = 2)

  axisLwd(1, line = -2)
  mtext(xlab, side = 1, line = 1.75, cex = par()$cex.lab * par()$cex)

}

#' Plot histogram
#'
#' Plot a histogram of the pre-industrial NGT-2012 values or trends and compare
#' it to the recent observations.
#'
#' @param piPeriod a vector of time points to specify the pre-industrial time
#'   period; defaults to the years 1000-1800 CE.
#' @param endRecentPeriod integer value for the last year to define the recent
#'   period for comparison with the pre-industrial period. The recent period is
#'   defined as the time period which ends in the year specified here and covers
#'   a total number of years corresponding to the size of the
#'   \code{filter.window}.
#' @param data.source character string to signal the data source: one of "iso"
#'   (d18O isotope data) or "acc" (accumulation data).
#' @param stack.method character; name of the NGT-2012 stacking method to use,
#'   see \code{selectHistogramData()} for details.
#' @param ... further parameters passed on to \code{selectHistogramData()} to
#'   control the stacking and merging methods.
#' @param filter.window single integer giving the size of the running mean
#'   window to use for filtering the data, and, if requested, the window over
#'   which linear trends are estimated; defaults to 11 (years).
#' @param type single character; must be either "anomaly" to plot the anomaly
#'   histograms, or "trend" to plot the linear trend histograms.
#' @param breaks a vector of break points between the histogram cells.
#' @param xlim x axis range; defaults to the range of \code{breaks}.
#' @param ylim y axis range.
#' @param col vector of length 2 specifying the colours to use to mark the
#'   pre-industrial and the recent data on the plot.
#' @param plot logical; set to \code{FALSE} to omit plotting of the histogram,
#'   so the function will only return the probability for the recent
#'   observations to occur under the pre-industrial distribution.
#' @param plot.quantiles logical; set to \code{TRUE} to plot vertical lines for
#'   the quantiles of the pre-industrial histogram specified by
#'   \code{quantile.probs}.
#' @param q.lty line type(s) of the vertical lines marking the quantiles set by
#'   \code{quantile.probs}; recycled to match the length of the latter.
#' @param quantile.probs the quantile probabilities of the pre-industrial
#'   histogram plotted when \code{plot.quantiles} is set to \code{TRUE}.
#' @param xmain optional main title above plot along x axis.
#' @param ymain optional main title next to y axis label.
#' @param plot.legend logical; shall a legend be plotted which marks the
#'   analysed time periods? Defaults to \code{TRUE}.
#' @param plot.2xaxes logical; shall a second x axis be plotted indicating an
#'   estimated temperature scale? Only used when isotope data is selected for
#'   the histogram analysis.
#' @param permil2temperature numeric; scale factor for the second x axis when
#'   \code{plot.2xaxes} is \code{TRUE}, i.e. the isotope-to-temperature slope.
#' @author Thomas Münch
#'
plotHistogram_DE <- function(piPeriod = 1000 : 1800, endRecentPeriod = 2011,
                             data.source = "iso", stack.method = "main", ...,
                             filter.window = 11, type = "anomaly",
                             breaks = seq(-2, 2, 0.2), xlim = range(breaks),
                             ylim = c(0, 1.25), col = c("black", "orange2"),
                             plot = TRUE, plot.quantiles = TRUE,
                             quantile.probs = c(0.005, 0.025, 0.975, 0.995),
                             q.lty = c(2, 4, 4, 2), xmain = NA, ymain = NA,
                             plot.legend = TRUE, plot.2xaxes = FALSE,
                             permil2temperature = 1 / 0.67) {

  if (!type %in% c("anomaly", "trend")) {
    stop("'type' must be either 'anomaly' or 'trend'.", call. = FALSE)
  }

  if (filter.window < 1) {
    stop("Running mean filter window needs to be >= 1.", call. = FALSE)
  }

  if ((filter.window %% 2) == 0) {
    stop("Running mean filter window needs to be odd.", call. = FALSE)
  }

  if (filter.window == 1 & type == "trend") {
    stop("Window size of 1 not suitable for trend estimation.", call. = FALSE)
  }

  # ----------------------------------------------------------------------------
  # Load data

  x <- selectHistogramData(data.source, stack.method, filter.window, ...)

  if (type == "anomaly") {

    var <- "stack"
    xlab <- grfxtools::LabelAxis(suffix = "anomaly")
    ylab <- grfxtools::LabelAxis("Relative Häufigkeitsdichte",
                                 unit.type = "freq")

    if (data.source == "acc") {

      xlab <- grfxtools::LabelAxis("Accumulation rate", unit = "mm w.eq.",
                                   unit.type = "trend", time.unit = "yr")
      ylab <- bquote("Probability density" * " (" * "mm w.eq."^{"-1"} ~ "yr)")
    }

    if (plot.2xaxes) {
      xlab <- "NGT-2012"
      pm.unit <- grfxtools::LabelAxis("")
      dc.unit <- grfxtools::LabelAxis("", unit = "celsius")
    }

  } else if (type == "trend") {

    x <- estimateSlopes(x, window = filter.window)

    var <- "slope"
    xlab <- grfxtools::LabelAxis(suffix = "trend", unit.type = "trend",
                                 time.unit = "yr")
    ylab <- grfxtools::LabelAxis("Probability density",
                                 unit = bquote("\u2030"^{"-1"} ~ "yr"))

    if (data.source == "acc") {

      xlab <- bquote("Accumulation rate trend" *
                     " (" * "mm w.eq." ~ "yr"^{"-2"} * ")")
      ylab <- bquote("Probability density" *
                     " (" * "mm w.eq."^{"-1"} ~ "yr"^{"2"} * ")")
    }

    if (plot.2xaxes) {
      xlab <- "NGT-2012 trend"
      pm.unit <- grfxtools::LabelAxis("", unit.type = "trend",
                                      time.unit = "yr")
      dc.unit <- grfxtools::LabelAxis("", unit.type = "trend", unit = "celsius",
                                      time.unit = "yr")
    }
  }

  # ----------------------------------------------------------------------------
  # Obtain distribution and recent data

  midpoint <- floor(filter.window / 2)
  t <- piPeriod[(midpoint + 1) : (length(piPeriod) - midpoint)]

  if (filter.window == 1) {
    recentTime <- seq(endRecentPeriod, by = -1, length.out = 11)
  } else {
    recentTime <- endRecentPeriod - floor(filter.window / 2)
  }

  piData <- subsetData(x, t, var)

  recentData <- x %>%
    subsetData(recentTime, var)

  # ----------------------------------------------------------------------------
  # Calculate probability for recent value to occur under pI distribution

  probability <- c(probability = pnorm(q = mean(recentData),
                                       mean = mean(piData), sd = sd(piData),
                                       lower.tail = FALSE))

  # ----------------------------------------------------------------------------
  # Make plot

  if (plot) {

    q <- seq(xlim[1], xlim[2], diff(breaks)[1] / 100)

    hst <- hist(piData, freq = FALSE, breaks = breaks, xlim = xlim, ylim = ylim,
                main = "", xlab = "", ylab = "", axes = FALSE,
                col = adjustcolor(col[1], alpha = 0.2), yaxs = "i")

    axx <- axisLwd(1)
    axy <- axisLwd(2)

    mtext(ylab, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0,
          adj = 0.215)

    if (plot.2xaxes) {

      at <- pretty(range(axx) * permil2temperature, length(axx))
      axis(1, line = 4.5, labels = FALSE, lwd = Lwd(1), lwd.ticks = 0)
      axis(1, at = at / permil2temperature, labels = at,
           line = 4.5, lwd = 0, lwd.ticks = Lwd(1))
      mtext(xlab, side = 1, line = 8.15, cex = par()$cex.lab * par()$cex)
      mtext(pm.unit, side = 1, line = 2.5,
            cex = par()$cex.lab * par()$cex, adj = 0.975)
      mtext(dc.unit, side = 1, line = 7,
            cex = par()$cex.lab * par()$cex, adj = 0.975)

    } else {

      mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)

    }

    lines(q, gq <- dnorm(q, mean = mean(piData), sd = sd(piData)),
          col = col[1], lwd = Lwd(2))

    ymax <- ifelse(max(hst$density) > max(gq), max(hst$density), max(gq))

    if (plot.quantiles) {
      q <- quantile(piData, quantile.probs, na.rm = TRUE)
      q.lty <- rep(q.lty, length.out = length(q))
      for (i in 1 : length(q)) {
        lines(rep(q[i], 2), c(0, ymax), col = col[1],
              lty = q.lty[i], lwd = Lwd(1))
      }
    }

    if (filter.window == 1) {

      par(new = TRUE)
      hist(recentData, freq = FALSE, breaks = breaks, xlim = xlim,
           axes = FALSE, main = "", xlab = "", ylab = "",
           col = adjustcolor(col[2], alpha = 0.2), yaxs = "i")

      ax <- axisLwd(4, col = col[2], col.axis = col[2])
      text(1.375 * xlim[2], y = mean(ax), "Probability density", srt = -90,
           xpd = NA, cex = par()$cex.lab * par()$cex, col = col[2])

    } else {

      lines(rep(recentData, 2), c(0, ymax), col = col[2], lwd = Lwd(3))

    }

    if (plot.legend) {

      recentString <- ifelse(filter.window == 1,
                             "Recent values", "Recent value")
      recentPeriod <- recentTime + (-midpoint : midpoint)

      lg <- c(sprintf("pI distribution\n(%s-%s CE)",
                      min(piPeriod), max(piPeriod)),
              "2.5, 50, 97.5 %\nquantiles",
              sprintf("%s\n(%s-%s CE)", recentString,
                      min(recentPeriod), max(recentPeriod)))
      lg.lty <- c(1, 5, 1)
      lg.lwd <- Lwd(c(10, 1, 2.5))
      lg.col <- c(adjustcolor(col[1], 0.2), col[1], col[2])

      if (type == "anomaly" & filter.window == 1) {
        lg.lwd[3] <- lg.lwd[1]
        lg.col[3] = adjustcolor(col[2], 0.2)
      }

      if (!plot.quantiles) {
        lg <- lg[-2]
        lg.lty <- lg.lty[-2]
        lg.lwd <- lg.lwd[-2]
        lg.col <- lg.col[-2]
      }

      legend("topleft", legend = lg, lty = lg.lty, lwd = lg.lwd, col = lg.col,
             bty = "n", adj = c(0, 0.75), y.intersp = 1.5, seg.len = 2.5)

    }

    if (!is.na(xmain)) {
      mtext(xmain, side = 3, line = 2.5, adj = 0,
            cex = 4 / 3 * par()$cex, font = 2)
    }
    if (!is.na(ymain)) {
      mtext(ymain, side = 2, line = 6, las = 0,
            cex = 4 / 3 * par()$cex, font = 2)
    }
  }

  invisible()

}

plot.NGT.MAR_DE <- function(filter.window = 11) {

  if (filter.window < 1) {
    stop("Running mean filter window needs to be >= 1.", call. = FALSE)
  }

  if ((filter.window %% 2) == 0) {
    warning("Running mean filter window should be odd.", call. = FALSE)
  }

  # ----------------------------------------------------------------------------
  # Load time series data

  permil2temperature.mid <- 1 / 0.67 # Greenland spatial slope
  permil2temperature.low <- 1 / 1.1  # Masson-Delmotte et al. (2015)
  permil2temperature.hig <- 2.1      # Vinther et al. (2009)

  filteredStackedNGT <- processNGT() %>%
    filterData(window = filter.window) %>%
    stackNGT() %>%
    dplyr::filter(Year >= 1871)

  filteredStackedTemperatureNGT.mid <- filteredStackedNGT %>%
    dplyr::mutate(stack = permil2temperature.mid * stack)
  filteredStackedTemperatureNGT.low <- filteredStackedNGT %>%
    dplyr::mutate(stack = permil2temperature.low * stack)
  filteredStackedTemperatureNGT.hig <- filteredStackedNGT %>%
    dplyr::mutate(stack = permil2temperature.hig * stack)

  filteredMAR <- readMAR() %>%
    filterData(window = filter.window)

  # ----------------------------------------------------------------------------
  # Make plots

  xlab  <- "Jahr n.d.Z."
  ylab1 <- "Anomalie des\nSchmelzwasserabfluss (Gt pro Jahr)"
  ylab2a <- "NGT-2012"
  ylab2b <- grfxtools::LabelAxis("Temperatur-Anomalie", unit = "celsius")

  xlim <- c(1870, 2012)
  ylim1 <- c(-75, 330)
  ylim2 <- c(-2.5, 3.5)

  x1 <- 2045
  x2 <- 2056
  y1 <- 0.5

  xanml <- c(1870.5,  2011.5)
  yanml <- c(0, 0)

  col <- c("black", "#d95f02")

  plot(filteredMAR$Year, filteredMAR$melt, type = "l", axes = FALSE,
       xlim = xlim, ylim = ylim1, xlab = "", ylab = "")

  lines(xanml, yanml, lty = 2, lwd = Lwd(1.5), col = col[2])
  lines(filteredMAR$Year, filteredMAR$melt, lwd = Lwd(3), col = col[2])

  axisLwd(1, at = seq(1870, 2010, 35), line = 0.5)
  axisLwd(2, at = seq(-50, 300, 50), col = col[2], col.axis = col[2])
  mtext(xlab, 1, 3.5, cex = par()$cex.lab)
  mtext(ylab1, 2, 3.6, col = col[2], cex = par()$cex.lab, las = 0,
        adj = 0.5)

  par(new = TRUE)

  plot(filteredStackedTemperatureNGT.mid, type = "n", axes = FALSE,
       xlab = "", ylab = "", xlim = xlim, ylim = ylim2)

  lines(xanml, yanml, lty = 2, lwd = Lwd(1.5), col = "darkgrey")
  rect(1870.5, -2.6, 2011.5, 3.25, border = "dimgrey", lwd = Lwd(2), xpd = NA)

  axisLwd(4)
  text(x2, y1, ylab2a, srt = -90, xpd = NA, cex = par()$cex.lab * par()$cex)
  text(x1, y1, ylab2b, srt = -90, xpd = NA, cex = par()$cex.lab * par()$cex)

  lines(filteredStackedTemperatureNGT.mid, lwd = Lwd(3), col = col[1])

  grfxtools::Polyplot(filteredStackedTemperatureNGT.mid$Year,
                      filteredStackedTemperatureNGT.low$stack,
                      filteredStackedTemperatureNGT.hig$stack,
                      col = col[1], alpha = 0.2)
  lines(filteredStackedTemperatureNGT.low, lwd = Lwd(1), col = col[1])
  lines(filteredStackedTemperatureNGT.hig, lwd = Lwd(1), col = col[1])

}

# ------------------------------------------------------------------------------
# Plot Reklim figure 2

asp <- 2 / 3
w <- 18.3
natfig(height = asp * w, width = w, file = "zzz/warnsignale_fig_test.pdf")

# sub-figure placement specifiers
x1 <- 0.5
x2 <- 0.6
y1 <- 0.3
y2 <- y1 - 0.05

# plot NGT time series
op <- par(mar = c(0, 0, 0, 0), oma = c(4.25, 5, 0, 5), fig = c(0, x1, y1, 1))
plotNGT()
mtext("a", side = 3, adj = -0.12, line = -1.5, font = 2, cex = 4 / 3)
mtext("b", side = 3, adj = 1.28, line = -1.5, font = 2, cex = 4 / 3)
par(op)

# plot histogram
op <- par(cex = 1, mar = c(9.5, 5, 1.5, 0.5), fig = c(x2, 1, y1, 1), new = TRUE)
plotHistogram_DE(type = "anomaly", plot.legend = FALSE, plot.2xaxes = TRUE,
                 q.lty = 2)
par(op)

# plot label for time series plot
op <- par(cex = 1, mar = c(0, 0, 0, 0), fig = c(0, x1, 0, y2), new = TRUE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

lg <- c("Vorindustrieller Zeitraum 1000\u20131800 n.d.Z.",
        "Rezenter Zeitraum 2001\u20132011 n.d.Z.")

legend("topleft", legend = lg, lty = c(1, 1), lwd = Lwd(c(20, 20)),
       col = c(adjustcolor("black", 0.2), adjustcolor("orange2", 0.8)),
       bty = "n", adj = c(0, 0.5), y.intersp = 2, seg.len = 2.5,
       inset = c(0.15, -0.075), cex = 7 / 6)
par(op)

# plot histogram label
op <- par(cex = 1, mar = c(0, 0, 0, 0), fig = c(x2, 1, 0, y2), new = TRUE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

lg <- c("Verteilung im vorindustriellen Zeitraum\n1000\u20131800 n.d.Z.\n",
        "Normalverteilung\n", "Quantile p = 0.95 und p = 0.99\n",
        "Mittelwert des rezenten Zeitraums\n2001\u20132011 n.d.Z.")

legend("topleft", legend = lg, lty = c(1, 1, 2, 1), lwd = Lwd(c(10, 2, 1, 3)),
       col = c(adjustcolor("black", 0.2), "black", "black", "orange2"),
       bty = "n", adj = c(0, 0.75), y.intersp = 0.5, seg.len = 2.5,
       inset = c(0.15, -0.15), cex = 7 / 6)
par(op)

dev.off()

# ------------------------------------------------------------------------------
# Plot Reklim figure 3a

w <- 8.5
h <- 5.25
natfig(height = h, width = w, file = "zzz/warnsignale_fig3a_test.pdf")

op <- par(mar = c(5, 7.25, 3, 7.25))
plot.NGT.MAR_DE(filter.window = 11)
mtext("c", side = 3, adj = -0.3, line = 2, font = 2, cex = 4 / 3)
par(op)

dev.off()
