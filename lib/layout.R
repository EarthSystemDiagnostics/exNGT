#' Produce figure plot in Nature layout
#'
#' This is a wrapper function for \code{grfxtools::Quartz} to produce plots in
#' the desired Nature format; including desired size in cm and desired font
#' sizes and line widths.
#'
#' @param file path to a file for storing a hardcopy of the plot including a
#'   supported file extension to set the \code{type} of output (e.g. ".pdf" or
#'   ".png"); see \code{?grfxtools::Quartz} for details.
#' @param height the height of the plotting area in cm.
#' @param width the width of the plotting area in cm.
#' @param ... further parameters passed on to \code{grfxtools::Quartz}.
#' @author Thomas Münch
natfig <- function(file = NULL, height, width, ...) {

  cm2inch <- 1 / 2.54
  h <- cm2inch * height
  w <- cm2inch * width

  grfxtools::Quartz(file = file, height = h, width = w,
                    cex.axis = 1, cex.lab = 7/6, cex.main = 7/6,
                    pointsize = 6, lwd = 0.44, ...)

}

#' Line width as a multiple of current \code{par} setting
#'
#' @param lwd numeric factor to scale current \code{par()$lwd} setting.
#' @return line width factor relative to the current \code{par()$lwd} setting.
#' @author Thomas Münch
Lwd <- function(lwd) {

  lwd * par()$lwd
}

#' Produce plot axis with line width controlled by current \code{par} setting
#'
#' @param ... parameters passed on to \code{axis()}.
#' @author Thomas Münch
axisLwd <- function(...) {

  axis(lwd = Lwd(1), ...)
}

#' Conversion from pt to mm using ggplot2 conversion factor
#'
#' @param pt size measured in pt.
#' @return size measured in mm for use in ggplot2 size specifications.
#' @author Thomas Münch
pt2mm <- function(pt) {

  pt / ggplot2::.pt

}
