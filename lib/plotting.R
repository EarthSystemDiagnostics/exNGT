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
plotHistogram <- function(piPeriod = 1000 : 1800, endRecentPeriod = 2011,
                          data.source = "iso", stack.method = "main", ...,
                          filter.window = 11, type = "anomaly",
                          breaks = seq(-2, 2, 0.2), xlim = range(breaks),
                          ylim = c(0, 1.25), col = c("black", "firebrick4"),
                          plot = TRUE, plot.quantiles = TRUE,
                          quantile.probs = c(0.025, 0.5, 0.975),
                          xmain = NA, ymain = NA, plot.legend = TRUE,
                          plot.2xaxes = FALSE, permil2temperature = 1 / 0.67) {

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
    ylab <- grfxtools::LabelAxis("Probability density", unit.type = "freq")

    if (data.source == "acc") {

      xlab <- grfxtools::LabelAxis("Accumulation rate", unit = "mm w.eq.",
                                   unit.type = "trend", time.unit = "yr")
      ylab <- bquote("Probability density" * " (" * "mm w.eq."^{"-1"} ~ "yr)")
    }

    if (plot.2xaxes) {
      xlab <- "NGT-2012 anomaly"
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

    axx <- axis(1)
    axy <- axis(2)

    mtext(ylab, side = 2, line = 3.25, cex = par()$cex.lab * par()$cex, las = 0)

    if (plot.2xaxes) {

      at <- pretty(range(axx) * permil2temperature, length(axx))
      axis(1, line = 5, labels = FALSE, lwd.ticks = 0)
      axis(1, at = at / permil2temperature, labels = at,
           line = 5, lwd = 0, lwd.ticks = 1)
      mtext(xlab, side = 1, line = 9.25, cex = par()$cex.lab * par()$cex)
      mtext(pm.unit, side = 1, line = 2.75,
            cex = par()$cex.lab * par()$cex, adj = 0.975)
      mtext(dc.unit, side = 1, line = 7.75,
            cex = par()$cex.lab * par()$cex, adj = 0.975)

    } else {

      axis(1)
      mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)

    }

    lines(q, gq <- dnorm(q, mean = mean(piData), sd = sd(piData)),
          col = col[1], lwd = 2)

    ymax <- ifelse(max(hst$density) > max(gq), max(hst$density), max(gq))

    if (plot.quantiles) {
      sapply(quantile(piData, quantile.probs, na.rm = TRUE), function(x) {
        lines(rep(x, 2), c(0, ymax), col = col[1], lty = 5, lwd = 1)
      })
    }

    if (filter.window == 1) {

      par(new = TRUE)
      hist(recentData, freq = FALSE, breaks = breaks, xlim = xlim,
           axes = FALSE, main = "", xlab = "", ylab = "",
           col = adjustcolor(col[2], alpha = 0.2), yaxs = "i")

      ax <- axis(4, col = col[2], col.axis = col[2])
      text(1.375 * xlim[2], y = mean(ax), "Probability density", srt = -90,
           xpd = NA, cex = par()$cex.lab * par()$cex, col = col[2])

    } else {

      lines(rep(recentData, 2), c(0, ymax), col = col[2], lwd = 2.5)

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
      lg.lwd <- c(10, 1, 2.5)
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
            cex = par()$cex.lab * par()$cex, font = 2)
    }
    if (!is.na(ymain)) {
      mtext(ymain, side = 2, line = 6, las = 0,
            cex = par()$cex.lab * par()$cex, font = 2)
    }
  }

  return(probability)

}

#' Plot map
#'
#' Plot a map of Greenland in polar projection with the NGT-2012 firn core and
#' Greenland weather station locations labelled.
#'
#' @author Thomas Münch
plotMap <- function() {

  cores <- loadPositions() %>%
    dplyr::filter(Identifier == "core")
  stations <- loadPositions() %>%
    dplyr::filter(Identifier == "station")

  min.lat <- 57.5
  max.lat <- 85
  min.lon <- -75
  max.lon <- -10

  lat.pos <- c(60, 70, 80)
  lon.pos <- -c(20, 40, 60)

  lat.pos.offset <- c(3.5, 5, 10)

  p <- grfxtools::ggpolar(pole = "N",
                          max.lat = max.lat, min.lat = min.lat,
                          max.lon = max.lon, min.lon = min.lon,
                          lat.ax.vals = lat.pos, long.ax.vals = lon.pos,
                          f.long.label.ticks = Inf, f.long.label.pos = 15,
                          rotate = TRUE, land.fill.colour = "transparent",
                          size.outer = 0.5,
                          lat.ax.labs.pos = min.lon - lat.pos.offset,
                          ax.labs.size = 4.75,
                          land.outline.colour = "burlywood4", clip = "off") +

  ggplot2::geom_text(data = cores, size = 2.5,
                     ggplot2::aes(x = Longitude, y = Latitude,
                                  label = Site)) +

  ggplot2::geom_label(data = stations,
                      ggplot2::aes(x = Longitude, y = Latitude,
                                   label = Site),
                      size = 2.5, alpha = 0.75, label.size = 0) +

  ggplot2::geom_point(data = cores,
                      ggplot2::aes(x = Longitude, y = Latitude),
                      col = "black", bg = "grey", size = 1.5,
                      pch = 21, stroke = 0.8) +

  ggplot2::geom_point(data = stations,
                      ggplot2::aes(x = Longitude, y = Latitude),
                      col = "black", size = 2.5, pch = 17)

  p

}

#' Produce a plot comparing NGT-2012 and Arctic2k time series
#'
#' @param filter.window single integer giving the size of the running mean
#'   window to use for filtering the isotope and Arctic2k data; defaults to 11
#'   (years).
#' @param permil2temperature numeric; isotope-to-temperature conversion slope
#'   (permil / K) to use for a second NGT-2012 plot axis; defaults to the
#'   Greenland spatial slope.
#' @author Thomas Münch
#'
plot.NGT.Arctic2k <- function(filter.window = 11,
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
    stackNGT()

  filteredStackedNGT <- NGT %>%
    filterData(window = filter.window) %>%
    stackNGT()

  Arctic2k <- readArctic2k() %>%
    extendWithHadCrut()

  filteredArctic2k <- Arctic2k %>%
    filterData(window = filter.window)

  # Linear regression data (annual means)

  t1 <- 1800 : 1000
  t2 <- 2011 : 1800

  regressionData <- list(
    data.frame(x = t1, y = subsetData(stackedNGT, t1, "stack")),
    data.frame(x = t2, y = subsetData(stackedNGT, t2, "stack")),
    data.frame(x = t1, y = subsetData(Arctic2k, t1, "TempAnomaly")),
    data.frame(x = t2, y = subsetData(Arctic2k, t2, "TempAnomaly"))
  )

  regressionModels <- regressionData %>%
    lapply(function(lst) {coef(lm(y ~ x, lst))})

  # ----------------------------------------------------------------------------
  # Make plots

  xlab  <- "Year CE"
  ylab.ngt.pm <- grfxtools::LabelAxis("NGT-2012")
  ylab.ngt.dc <- grfxtools::LabelAxis("NGT-2012", unit = "celsius")
  ylab.a2k <- grfxtools::LabelAxis("Arctic 2k", unit = "celsius")

  xlim <- c(1000, 2020)
  ylim.ngt <- c(-6, 2.5)
  ylim.a2k <- c(-2, 10.75)

  xanml <- c(800, 2020)
  yanml <- rep(0, 2)

  startNew <- 1993
  i <- match(startNew, stackedNGT$Year)
  n <- nrow(stackedNGT)
  startHadCrut <- 2000
  j <- match(startHadCrut, Arctic2k$Year)
  m <- nrow(Arctic2k)

  x1 <- 845
  x2 <- 2185
  y1 <- 0.
  y2 <- 0.5

  col <- c("black", "dodgerblue4")

  op <- par(mar = c(0, 0, 0, 0), oma = c(5, 5, 0, 5))

  plot(stackedNGT, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = xlim, ylim = ylim.ngt)

  lines(xanml, yanml, lty = 2, lwd = 1.5, col = "darkgrey")

  lines(stackedNGT, col = "darkgrey")

  lines(filteredStackedNGT[i : n, ], col = col[1], lwd = 2.5)
  lines(filteredStackedNGT[1 : i, ], col = "firebrick3", lwd = 2.5)

  axis(2, at = seq(-2, 2, 1))
  axis(4, labels = seq(-2, 3, 1), at = seq(-2, 3, 1) / permil2temperature)

  text(x1, y1, ylab.ngt.pm, srt = +90, xpd = NA,
       cex = par()$cex.lab * par()$cex, col = col[1])
  text(x2, y2 / permil2temperature, ylab.ngt.dc, srt = -90, xpd = NA,
       cex = par()$cex.lab * par()$cex, col = col[1])

  mtext("a", side = 3, adj = 0.01, line = -2.5, font = 2, cex = par()$cex.lab)

  lines(t1, regressionModels[[1]][1] + regressionModels[[1]][2] * t1,
        col = col[1], lwd = 2, lty = 2)
  lines(t2, regressionModels[[2]][1] + regressionModels[[2]][2] * t2,
        col = col[1], lwd = 2, lty = 2)

  par(new = TRUE)

  plot(Arctic2k$Year, Arctic2k$TempAnomaly, type = "n", axes = FALSE,
       xlab = "", ylab = "", xlim = xlim, ylim = ylim.a2k)

  axis(1)
  axis(4, at = seq(-2, 3, 1), col = col[2], col.axis = col[2])

  mtext(xlab, side = 1, line = 3.5, cex = par()$cex.lab * par()$cex)
  text(x2, y2, ylab.a2k, srt = -90, xpd = NA,
       cex = par()$cex.lab * par()$cex, col = col[2])

  lines(xanml, yanml, lty = 2, lwd = 1.5, col = "darkgrey")

  lines(Arctic2k$Year, Arctic2k$TempAnomaly,
        col = adjustcolor(col[2], alpha = 0.6))

  lines(filteredArctic2k$Year[1 : j], filteredArctic2k$TempAnomaly[1 : j],
        col = "deepskyblue1", lwd = 2.5)
  lines(filteredArctic2k$Year[j : m], filteredArctic2k$TempAnomaly[j : m],
        col = col[2], lwd = 2.5)

  lines(t1, regressionModels[[3]][1] + regressionModels[[3]][2] * t1,
        col = col[2], lwd = 2, lty = 2)
  lines(t2, regressionModels[[4]][1] + regressionModels[[4]][2] * t2,
        col = col[2], lwd = 2, lty = 2)

  mtext("b", side = 3, adj = 0.01, line = -20, font = 2, cex = par()$cex.lab)

  par(op)

}

#' Produce a plot comparing NGT-2012 and MAR time series
#'
#' @param filter.window single integer giving the size of the running mean
#'   window to use for filtering the isotope and Arctic2k data; defaults to 11
#'   (years).
#' @author Thomas Münch
#'
plot.NGT.MAR <- function(filter.window = 11) {

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
    stackNGT()

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

  xlab  <- "Year CE"
  ylab1 <- grfxtools::LabelAxis("Temperature", unit = "celsius")
  ylab2 <- grfxtools::LabelAxis("Melt runoff", unit = "Gt",
                                unit.type = "trend", time.unit = "yr")

  xlim <- c(1875, 2011)
  ylim1 <- c(-2.5, 3.5)
  ylim2 <- c(-55, 300)

  x1 <- 2034
  y1 <- 125

  col <- c("black", "#d95f02", "#7570b3")

  plot(filteredStackedTemperatureNGT.mid, type = "n", axes = FALSE,
       xlab = "", ylab = "", xlim = xlim, ylim = ylim1)

  axis(1, at = seq(1870, 2010, 20))
  axis(2)
  mtext(xlab, 1, 3.5, cex = par()$cex.lab)
  mtext(ylab1, 2, 3, cex = par()$cex.lab, las = 0)
  mtext("c", side = 3, adj = 0.01, line = -1, font = 2, cex = par()$cex.lab)

  lines(filteredStackedTemperatureNGT.mid, lwd = 2, col = col[1])

  grfxtools::Polyplot(filteredStackedTemperatureNGT.mid$Year,
                      filteredStackedTemperatureNGT.low$stack,
                      filteredStackedTemperatureNGT.hig$stack,
                      col = col[1], alpha = 0.2)
  lines(filteredStackedTemperatureNGT.low, lwd = 1, col = col[1])
  lines(filteredStackedTemperatureNGT.hig, lwd = 1, col = col[1])

  lines(filteredMAR$Year, filteredMAR$t2m, lwd = 2, col = col[2])

  par(new = TRUE)

  plot(filteredMAR$Year, filteredMAR$melt, type = "l", axes = FALSE,
       col = col[3], lwd = 2, xlim = xlim, ylim = ylim2,
       xlab = "", ylab = "")

  axis(4, col = col[3], col.axis = col[3])
  text(x1, y1, ylab2, srt = -90, xpd = NA,
       col = col[3], cex = par()$cex.lab * par()$cex)

  legend("topleft",
         c("NGT-2012 temperature", "MAR3.5 temperature", "MAR3.5 melt runoff"),
         inset = c(0, 0.05),
         lty = 1, lwd = 2, col = col, bty = "n")

}

#' Produce spectrum plot of NGT-2012 versus Arctic2k
#'
#' @param filter.window single integer giving the size of the running mean
#'   window to use for filtering the isotope and Arctic2k data; defaults to 11
#'   (years).
#' @author Thomas Münch
#'
plotSpectrum <- function(filter.window = 11) {

  if (filter.window < 1) {
    stop("Running mean filter window needs to be >= 1.", call. = FALSE)
  }

  if ((filter.window %% 2) == 0) {
    warning("Running mean filter window should be odd.", call. = FALSE)
  }

  # ----------------------------------------------------------------------------
  # Load data

  # Parameter - convert NGT to temperature with three different calibrations

  permil2temperature.mid <- 1 / 0.67 # Greenland spatial slope
  permil2temperature.low <- 1 / 1.1  # Masson-Delmotte et al. (2015)
  permil2temperature.hig <- 2.1      # Vinther et al. (2009)

  # Spectra

  spectrumA2k <- readArctic2k() %>%
    subsetData(t = 1979 : 1505, var = "TempAnomaly") %>%
    ts(deltat = 1) %>%
    proxysnr:::SpecMTM() %>%
    proxysnr:::LogSmooth()

  spectraNGT <- list()

  spectraNGT$mid <- selectNGTForSpectra() %>%
    proxysnr::ArraySpectra(df.log = 0.1) %>%
    proxysnr::SeparateSpectra() %>%
    .$signal

  spectraNGT$low <- spectraNGT$hig <- spectraNGT$mid

  # apply calibrations
  spectraNGT$mid$spec <- spectraNGT$mid$spec * (permil2temperature.mid^2)
  spectraNGT$low$spec <- spectraNGT$low$spec * (permil2temperature.low^2)
  spectraNGT$hig$spec <- spectraNGT$hig$spec * (permil2temperature.hig^2)

  # ----------------------------------------------------------------------------
  # Plot

  col <- c("black", "dodgerblue4")

  xlab  <- "Time period (yr)"
  ylab <- grfxtools::LabelAxis("Power spectral density", unit = "celsius",
                               time.unit = "yr", unit.type = "psd")
  tcklab <- c(expression(10^{-1}), expression(10^{0}), expression(10^{1}))

  n.crit.lower <- 1

  proxysnr:::LPlot(spectraNGT$mid, bPeriod = TRUE, bNoPlot = TRUE, axes = FALSE,
                   xlab = "", ylab = "", xlim = c(225, 5), ylim = c(0.05, 10),
                   xaxs = "i")

  axis(1)
  axis(2, at = c(0.1, 1, 10), labels = tcklab)

  mtext(xlab, 1, 3.5, cex = par()$cex.lab)
  mtext(ylab, 2, 3.5, cex = par()$cex.lab, las = 0, adj = 0.925)

  n.crit.upper <- length(which(spectraNGT$mid$freq > 1 / 5))
  proxysnr:::LLines(spectraNGT$mid, bPeriod = TRUE, lwd = 3, col = col[1],
                    removeFirst = n.crit.lower, removeLast = n.crit.upper)

  i.keep <- (n.crit.lower + 1) : (length(spectraNGT$mid$freq) - n.crit.upper)
  proxysnr:::Polyplot(1 / spectraNGT$mid$freq[i.keep],
                      spectraNGT$hig$spec[i.keep], spectraNGT$low$spec[i.keep],
                      col = col[1], alpha = 0.2)
  proxysnr:::LLines(spectraNGT$low, bPeriod = TRUE, lwd = 1, lty = 1,
                    removeFirst = n.crit.lower, removeLast = n.crit.upper,
                    col = col[1])
  proxysnr:::LLines(spectraNGT$hig, bPeriod = TRUE, lwd = 1, lty = 1,
                    removeFirst = n.crit.lower, removeLast = n.crit.upper,
                    col = col[1])

  n.crit.upper <- length(which(spectrumA2k$freq > 1 / 5))
  proxysnr:::LLines(spectrumA2k, bPeriod = TRUE, conf = FALSE, lwd = 3,
                    removeFirst = n.crit.lower, removeLast = n.crit.upper,
                    col = col[2])

  legend("bottomleft", c("NGT-2012", "Arctic 2k"),
         lty = 1, lwd = 3, col = col, bty = "n")

  # ----------------------------------------------------------------------------
  # How much higher is NGT variability?

  combinedSpectra <- data.frame(
    freq = spectrumA2k$freq,
    a2k = spectrumA2k$spec,
    ngtMid = spectraNGT$mid$spec,
    ngtHig = spectraNGT$hig$spec,
    ngtLow = spectraNGT$low$spec)

  timescale <- c(1 / 51, 1 / 11)

  variabilityRatio <- combinedSpectra %>%
    dplyr::filter(freq >= timescale[1] & freq <= timescale[2]) %>%
    dplyr::summarise(dplyr::across(.fns = mean)) %>%
    dplyr::transmute(ngtLow = ngtLow / a2k,
                     ngtMid = ngtMid / a2k,
                     ngtHig = ngtHig / a2k)

  cat("\nVariability ratio NGT vs. A2k (11 - 51 yr time scale):\n")
  print(variabilityRatio)
  cat("\n")

}

# ------------------------------------------------------------------------------
# main figures

#' Produce paper figure 01
#'
#' @param filter.window single integer giving the size of the running mean
#'   window to use for filtering the isotope and Arctic2k data; defaults to 11
#'   (years).
#' @param permil2temperature numeric; isotope-to-temperature conversion slope
#'   (permil / K) to use for a second NGT-2012 plot axis; defaults to the
#'   Greenland spatial slope.
#' @author Thomas Münch
#'
makeFigure01 <- function(filter.window = 11, permil2temperature = 1 / 0.67) {

  y1 <- 0.41
  y2 <- 0.35

  par(fig = c(0, 1, y1, 1))
  plot.NGT.Arctic2k(filter.window = filter.window,
                    permil2temperature = permil2temperature)

  par(mar = c(5, 5, 0.5, 5), fig = c(0, 1, 0, y2), new = TRUE)
  plot.NGT.MAR(filter.window = filter.window)
  lines(x = c(2011, 2011), y = c(425, 280), col = "darkgrey",
        lwd = 1.5, xpd = NA)
  ## lines(x = c(1991, 1991), y = c(425, 375), col = "darkgrey",
  ##       lwd = 1.5, xpd = NA)
  ## lines(x = c(1991, 1871), y = c(375, 280), col = "darkgrey",
  ##       lwd = 1.5, xpd = NA)
  lines(x = c(1871, 1941), y = c(280, 364.5833), col = "darkgrey",
        lwd = 1.5, xpd = NA)
  lines(x = c(1956, 1991), y = c(382.7083, 425), col = "darkgrey",
        lwd = 1.5, xpd = NA)

}

#' Produce paper figure 02
#'
#' @param filter.window single integer giving the size of the running mean
#'   window to use for filtering the isotope and Arctic2k data; defaults to 11
#'   (years).
#' @author Thomas Münch
#'
makeFigure02 <- function(filter.window = 11) {

  plotSpectrum(filter.window = filter.window)

}

#' Produce paper figure 03
#'
#' @author Thomas Münch
#'
makeFigure03 <- function() {

  layout(matrix(1 : 2, 1, 2), widths = c(0.7, 0.3))
  par(cex = 1, mar = c(11, 5, 0.5, 0.5))

  plotHistogram(type = "anomaly", plot.legend = FALSE, plot.2xaxes = TRUE)

  par(mar = c(0, 0, 0, 0))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  lg <- c("Pre-ind. distribution\n(1000-1800 CE)", "Gaussian fit",
          "2.5, 50, 97.5 %\nquantiles", "Recent value\n(2001-2011 CE)")

  legend("topleft", legend = lg, lty = c(1, 1, 5, 1), lwd = c(10, 2, 1, 2.5),
         col = c(adjustcolor("black", 0.2), "black", "black", "firebrick4"),
         bty = "n", adj = c(0, 0.75), y.intersp = 1.5, seg.len = 2.5)

}
