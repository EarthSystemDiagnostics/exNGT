#' Calculate the diffusivity in polar firn.
#'
#' This function calculates the diffusivity in polar firn for the stable water
#' isotopes oxygen-18 and deuterium, depending on site-specific parameters.
#'
#' The calculations are primarily based on the following expressions in Johnsen
#' et al. (2000): The general expression for the diffusivity in firn of
#' oxygen-18 and deuterium is given by Eq. (17). The tortuosity of the firn is
#' accounted for by applying Eq. (18) with a tortuosity constant b = 1.3.  For
#' the diffusivity of water vapour in air we use Eq. (19). The temperature
#' dependence of the fractionation factors as well as the diffusivities in air
#' for oxygen-18 and deuterium are given in the text on p. 127. Note that here,
#' Johnsen et al. accidentally mix the factors controlling the isotope
#' diffusivities in air (compare the original reference Merlivat and Jouzel,
#' 1979). This has been accounted for. For the water vapour saturation vapour
#' pressure over ice we use, instead of the expression given on p. 126, the
#' parameterization given in van der Wel et al. (2015) (Eq. 5).
#'
#' @references
#' Johnsen, S. J., Clausen, H. B., Cuffey, K. M., Hoffmann, G., Schwander, J.,
#' and Creyts, T.: Diffusion of stable isotopes in polar firn and ice: the
#' isotope effect in firn diffusion, in: Physics of ice core records, edited
#' by: Hondoh, T., vol. 159, Hokkaido Univ. Press, Sapporo, Japan, 121–140,
#' 2000.
#' 
#' Merlivat, L. and J. Jouzel: Global climatic interpretation of the
#' deuterium-oxygen 18 relationship for precipitation. J. Geophys. Res.:
#' Oceans, 84, C8, 5029-5033, doi: 10.1029/JC084iC08p05029, 1979.
#' 
#' van der Wel, L. G., H. A. Been, R. S. W. van de Wal, C. J. P. P. Smeets and
#' H. A. J. Meijer: Constraints on the d2H diffusion rate in firn from
#' field measurements at Summit, Greenland. The Cryosphere, 9, 1089–1103, doi:
#' 10.5194/tc-9-1089-2015, 2015.
#'
#' @param rho Numeric vector of firn densities [kg/m^3] for which diffusivity
#'   is calculated. 
#' @param T firn temperature in [K].
#' @param P local surface pressure in [mbar].
#' @param dD if \code{TRUE} the diffusivity for deuterium is returned,
#'   otherwise for oxygen-18. Defaults to \code{FALSE}.
#' @return A numeric vector of calculated diffusivities in [cm^2/s] for the
#'   firn densities given by \code{rho}.
#' @author Thomas Muench, modified by Thomas Laepple
Diffusivity <- function(rho, T, P, dD = FALSE) {
  
  # Set physical constants
  kR <- 8.314478               # Gas constant [J/(K * mol)]
  kM <- 18.02e-3               # molar weight of H2O molecule [kg/mol]
  kP0 <- 1013.25               # standard atmospheric pressure [mbar]
  kRhoIce <- 920.              # density of ice [kg/m3]

  # Saturation vapour pressure over ice [Pa]
  p <- exp(9.5504 + 3.53 * log(T) - 5723.265 / T - 0.0073 * T)
  # Tortuosity constant
  b <- 1.3

  # Set fractionation factor
  if (dD) {
    alpha <- exp(16288 / (T^2) - 9.45e-2)
  } else {
    alpha <- exp(11.839 / T - 28.224e-3)
  }
  
  # Calculate water vapour diffusivity in air [m^2/s]
  Da <- 2.11e-5 * (T / 273.15)^(1.94) * (kP0 / P)
  if (dD) {
    isotopeFactor <- 1.0251
  } else {
    isotopeFactor <- 1.0285
  }
  Dai <- Da / isotopeFactor

  # Calculate tortuosity
  invtau <- rep(NA, length(rho))
  for (i in 1 : length(rho))
  {
    if (rho[i] <= kRhoIce / sqrt(b)) {
      invtau[i] <- 1 - b * (rho[i] / kRhoIce)^2
    } else {
      invtau[i] <- 0
    }
  }

  # Calculate isotope diffusivity in firn [m^2/s]
  D <- (kM * p * invtau * Dai * (1 / rho - 1 / kRhoIce)) / (kR * T * alpha)

  # Return firn diffusivity in [cm^2/s]
  return(D * 1e4)
  
}

#' Analytical solution of Herron-Langway firn densification model.
#'
#' This function calculates firn density depending on depth based on the
#' analytical solution of the steady-state Herron-Langway densification model.
#'
#' The empirical, steady-state Herron-Langway [HL] model of firn densification
#' together with its analytical solution is described in Herron and Langway
#' (1980). Its implementation here is based on the MATLAB code by Aslak
#' Grinsted (2014) which follows the nomenclature by Arthern et al. (2010). For
#' details of expressing the analytical solution, see also van der Wel (2012).
#'
#' The HL model is matched to density observations from Greenlandic and
#' Antarctic firn cores. The correction factors of Johnsen et al. (2000) (set
#' for \code{JohnsenCorr = TRUE}) have been introduced to further improve the
#' match with central Greenland firn core data.
#'
#' @references
#' Herron, M. M. and Langway Jr., C. C.: Firn densification: an empirical
#' model, J. Glaciol., 25(93), 373-385, 1980.
#'
#' Grinsted, A.: Steady state snow and firn density model,
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/47386-steady-state-snow-and-firn-density-model/content//densitymodel.m}, 2014. 
#'
#' Arthern, J. A., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas,
#' E. R.: In situ measurements of Antarctic snow compaction compared with
#' predictions of models, J. Geophys. Res., 115(F03011), 2010.
#'
#' van der Wel, L. G.: Analyses of water isotope diffusion in firn:
#' contributions to a better palaoclimatic interpretation of ice cores, Doctor
#' of Philosophy, University of Groningen,
#' \url{http://www.rug.nl/research/portal/files/2408704/Thesis.pdf}, 2012.
#' 
#' Johnsen, S. J., Clausen, H. B., Cuffey, K. M., Hoffmann, G., Schwander, J.,
#' and Creyts, T.: Diffusion of stable isotopes in polar firn and ice: the
#' isotope effect in firn diffusion, in: Physics of ice core records, edited
#' by: Hondoh, T., vol. 159, Hokkaido Univ. Press, Sapporo, Japan, 121–140,
#' 2000.
#'
#' @param depth Numeric vector of firn depths [m] at which firn density is
#'     calculated.
#' @param rho.surface surface density in [kg/m^3].
#' @param T mean firn temperature in [K].
#' @param bdot local mass accumulation rate in [kg/m^2/year].
#' @param JohnsenCorr logical; whether or not to apply the Johnsen correction
#'   to the Arrhenius rate constants for central Greenland sites.
#' @return A list with two elements:
#'     \itemize{
#'     \item \code{depth.we}: Numeric vector of water-equivalent depth in [m]
#'     corresponding to the true firn depths given by \code{depth}.
#'     \item \code{rho}: Numeric vector (length of \code{depth}) of firn
#'     density in [kg/m^3].}
#' @author Thomas Laepple, modified by Thomas Muench
DensityHL <- function(depth = (0 : 9000) / 100, rho.surface, T, bdot,
                      JohnsenCorr = FALSE) {

  # Constants
  kRhoIce <- 920.      # Density of ice [kg/m^3]
  kRhoW <- 1000.       # Density of water [kg/m^3]
  kR <- 8.314478       # Gas constant [J/(K * mol)]

  # Critical density of point between settling and creep-dominated stages
  kRhoC <- 550

  # Herron-Langway Arrhenius rate constants
  k0 <- 11 * exp(-10160 / (kR * T))
  k1 <- 575 * exp(-21400 / (kR * T))

  # Johnsen et al. (2000) correction for central Greenland cores
  if (JohnsenCorr) {
    k0 <- 0.85 * k0
    k1 <- 1.15 * k1
  }

  # Rate constants for time-dependent densification
  # (original Eq. (4) in Herron and Langway et al. (1980))
  A <- bdot / kRhoW
  c0 <- k0 * A
  c1 <- k1 * sqrt(A)

  # Rate constants for depth-dependent steady-state densification
  # (from converting the full time derivative to a depth derivative
  # neglecting the partial time derivative to get steady-state solution)
  d0 <- c0 / bdot
  d1 <- c1 / bdot

  fac.r0 <- rho.surface / (kRhoIce - rho.surface)
  fac.rc <- kRhoC / (kRhoIce - kRhoC)
  
  # Critical depth at which density reaches kRhoC
  z.c <- log(fac.rc / fac.r0) / (kRhoIce * d0)

  index.upper <- which(depth <= z.c)
  index.lower <- which(depth > z.c)

  # Steady-state density profile
  q <- rep(NA, length(depth))
  q[index.upper] <- fac.r0 * exp(d0 * kRhoIce * depth[index.upper])
  q[index.lower] <- fac.rc * exp(d1 * kRhoIce* (depth[index.lower] - z.c))
  rho <- kRhoIce * (q / (1+q))

  # Time when critical depth is reached
  tmp <- (kRhoIce - rho.surface) / (kRhoIce - kRhoC)
  t.c <- log(tmp) / c0
  
  # Steady-state time - water-equivalent depth relation
  t <- rep(NA, length(depth))
  tmp <- (kRhoIce - rho.surface) / (kRhoIce - rho[index.upper])
  t[index.upper] <- log(tmp) / c0
  tmp <- (kRhoIce - kRhoC) / (kRhoIce - rho[index.lower])
  t[index.lower] <- log(tmp) / c1 + t.c

  depth.we <- A * t

  return(list(depth.we = depth.we, rho = rho))

}

#' Calculate the diffusion length in polar firn.
#'
#' This function calculates the diffusion length in polar firn for the stable
#' water isotopes oxygen-18 and deuterium, depending on site-specific
#' parameters.
#'
#' The calculation of the diffusion length in firn is an implementation of
#' Eq. (8) in Gkinis et al. (2014) and is partly inspired by the corresponding
#' implementation of the PRYSM model
#' (\url{https://github.com/sylvia-dee/PRYSM}) presented in Dee et al. (2015).
#'
#' As input, a depth and a density vector have to be provided. If only a single
#' density value is passed to the function, the function silently repeats this
#' density value to build a density vector that matches the length of
#' \code{depth} and calculates the diffusion length for zero strain rate (not
#' yet implemented!). For a single temperature value as input, one diffusivity
#' value from calling the \code{\link{Diffusivity}} function is used to
#' calculate the diffusion length. Providing a vector of temperatures (which
#' length has to match \code{depth}, otherwise the function exits with an
#' error), results in polythermal diffusivity calculation where for each set of
#' depth, density and temperature, the diffusivity is calculated. This
#' depth-dependent diffusivity is then used to calculate the diffusion
#' lengths.
#'
#' \code{bFill} controls the output of the final diffusion length value at the
#' bottom of \code{depth}. This value depends on the unknown density and
#' related gradients at this position. For \code{bFill = TRUE} (the default),
#' the last known gradients are used for the calculation of the final diffusion
#' length value. Otherwise \code{NA} is returned.
#'
#' @references
#' Gkinis, V., Simonsen, S. B., Buchardt, S. L., White, J. W. C., and Vinther,
#' B. M.: Water isotope diffusion rates from the North-GRIP ice core for the
#' last 16,000 years – Glaciological and paleoclimatic implications, Earth
#' Planet. Sc. Lett., 405, 132–141, doi:10.1016/j.epsl.2014.08.022, 2014.
#' 
#' Dee, S., Emile-Geay, J., Evans, M. N., Allam, A., Steig, E. J., and
#' Thompson, D. M.: PRYSM: An open-source framework for PRoxY System Modeling,
#' with applications to oxygen-isotope systems, J. Adv. Model. Earth Syst., 7,
#' 1220–1247, doi:10.1002/2015MS000447, 2015.
#'
#' @param depth Numeric vector of firn depths [m] at which the diffusion
#'   lengths are calculated.
#' @param rho Numeric vector of firn density [kg/m^3], either of length one or
#'   of same length as \code{depth}.
#' @param T Numeric vector of firn temperature [K], either of length one or of
#'   same length as \code{depth}. Defaults to 10 m firn temperature at Kohnen
#'   Station.
#' @param P local surface pressure in [mbar]. Defaults to mean AWS9 value at
#'   Kohnen Station.
#' @param bdot local mass accumulation rate in [kg/m^2/year]. Defaults to
#'   long-time mean value at Kohnen Station.
#' @param dD if \code{TRUE} the diffusion length for deuterium is returned,
#'   otherwise for oxygen-18. Defaults to \code{FALSE}.
#' @param bFill if \code{TRUE} (the default) use the last known density
#'   and related gradients for the value at the bottom of \code{depth} to
#'   calculate the final diffusion length; see Details.
#' @return Numeric vector of the calculated diffusion lengths in [cm] at the
#'   depths given by \code{depth}.
#' @author Thomas Muench, modified by Thomas Laepple
DiffusionLength <- function(depth, rho, T = 273.15 - 44.5, P = 677, bdot = 64,
                            dD = FALSE, bFill = TRUE) {
  
  # Density of water
  kRhoW <- 1000.

  z <- depth

  # TODO (tmuench): implement zero-strain rate solution
  #if (length(rho) == 1) rho <- rep(rho, length(z))
  if (length(rho) != length(z))
    stop("Conflicting INPUT: 'depth' and 'rho' have different lengths")
  
  # Depth increments
  # CHECK (tlaepple): get dz in (m) from the z vector (in m) + extend with the
  # mean
  dz <- c(diff(z), mean(diff(z)))
  
  # Set time scale accounting for densification
  time_d <- cumsum(dz / (bdot / kRhoW) * (rho / kRhoW))
  # Convert from years to seconds
  ts <- time_d * 365.25 * 24 * 3600

  # Approximate density and related gradients
  drho <- diff(rho)
  dtdrho <- diff(ts) / diff(rho)

  # Fill unknown gradients at final depth
  ifelse(bFill,
         fill.gradient <- c(drho[length(drho)], dtdrho[length(dtdrho)]),
         fill.gradient <- rep(NA, 2))
  drho <- c(drho, fill.gradient[1])
  dtdrho <- c(dtdrho, fill.gradient[2])

  # Calculate diffusivity
  D <- vector(mode = "numeric", length = length(rho))
  if (length(T) == 1) {
    D <- Diffusivity(rho, T, P, dD = dD)
  } else {
    if (length(T) != length(rho))
      stop(paste("T and depth must have the same length",
                 "for polythermal diffusivity calculation"))
    for (i in 1 : length(rho))
      D[i] <- Diffusivity(rho[i], T[i], P, dD = dD)
  }
  
  # Integrate diffusivity along the density gradient
  # to obtain diffusion length [cm]
  sigma_sqrd_dummy <- 2 * (rho^2) * dtdrho * D
  sigma_sqrd <- cumsum(sigma_sqrd_dummy * drho)
  sigma <- sqrt(1 / (rho^2) * sigma_sqrd)

  return(sigma)

}

#' Calculate the diffusion length in years
#'
#' This function calculates the diffusion length in firn in time units (years)
#' based on calculating the diffusion length in depth units
#' (see \code{\link{DiffusionLength}}) and converting it from depth units in
#' temporal units using the firn density from the Herron-Langway model
#' (see \code{\link{DensityHL}}).
#'
#' The diffusion length can be calculated for several sites with varying
#' climatic input parameters \code{T}, \code{P}, \code{bdot} and
#' \code{rho.surface}. Note that for this all input parameter vectors must
#' either have the same length. Else, length-one vectors are recycled to match
#' the length of the longest input; if this still results in varying vector
#' lengths, an error is issued.
#'
#' @param core.length the simulated core length in metre for calculating the
#'   Herron-Langway firn density and the diffusion length. If this length is
#'   not sufficient to cover the requested time span given by \code{nt} and
#'   \code{t.res}, diffusion length values for the remaining time points are
#'   filled with the last properly obtained value. This issues only a warning
#'   since if the simulated core is long enough to reach the ice, the diffusion
#'   length is anyway constant, but it is a problem if the simulated core is
#'   far too shallow.
#' @param z.res the resolution in metre of the simulated firn core.
#' @param nt the number of time points for the output temporal diffusion
#'   length.
#' @param t.res the temporal resolution for the diffusion length in time units
#'   [yr]; i.e. the total time span covered is \code{t.res * nt}.
#' @param T local annual mean surface temperature (10-m firn temperature) in
#'   [K].
#' @param P local surface pressure in [mbar].
#' @param bdot local mass accumulation rate in [kg/m^2/yr].
#' @param rho.surface local surface firn density in [kg/m^3].
#' @param dD logical; if \code{TRUE} the diffusion length for deuterium is
#'   returned, otherwise for oxygen-18. Defaults to \code{FALSE}.
#' @param JohnsenCorr logical; whether or not to apply the Johnsen correction
#'   to the Arrhenius rate constants in the Herron-Langway model for central
#'   Greenland sites.
#' @param names optional character vector of site names.
#' @return the temporal diffusion lengths in a data frame of \code{nt} rows and
#'   a minimum of two columns, where the first column is the time axis and the
#'   second or other columns are the diffusion lengths for the requested
#'   site(s).
#' @author Thomas Münch
TemporalDiffusionLength <- function(core.length = 1000, z.res = 0.01,
                                    nt = 2000, t.res = 1,
                                    T = 273.15 - 44.5, P = 677,
                                    bdot = 64, rho.surface = 340,
                                    dD = FALSE, JohnsenCorr = FALSE,
                                    names = NULL) {

  nn <- c(length(T), length(P), length(bdot), length(rho.surface))
  if (stats::sd(nn) > 0) {

    if (nn[1] == 1) T <- rep(T, max(nn))
    if (nn[2] == 1) P <- rep(P, max(nn))
    if (nn[3] == 1) bdot <- rep(bdot, max(nn))
    if (nn[4] == 1) rho.surface <- rep(rho.surface, max(nn))

    nn <- c(length(T), length(P), length(bdot), length(rho.surface))

    if (stats::sd(nn) > 0) {
      stop("Inconsistent input length of 'T', 'P', 'bdot', 'rho.surface'.",
           call. = FALSE)
    }
  }

  n <- nn[1]

  # Density of water
  kRhoW <- 1000.

  # Depth profile
  depth <- seq(z.res, core.length, z.res)

  # Equidistant age scale
  t.equi <- seq_len(nt) * t.res

  # Loop over individual core sites
  sigma <- sapply(1 : n, function(i) {

    # Herron-Langway density profile for site 'i'
    HL  <- DensityHL(depth = depth,
                     rho.surface = rho.surface[i], T = T[i], bdot = bdot[i],
                     JohnsenCorr = JohnsenCorr)

    # Age of core according to Herron-Langway solution and const. acc. rate
    t <- HL$depth.we * kRhoW / bdot[i]

    # Check simulated vs. requested age
    if (max(t) < nt * t.res) {
      warning("Max. of simulated age too small; increase core length.",
              call. = FALSE)
    }

    # Diffusion length in [cm] for site 'i' as a function of depth
    sig.cm <- DiffusionLength(depth = depth, rho = HL$rho, T = T[i],
                              P = P[i], bdot = bdot[i], dD = dD)

    # convert diffusion length from [cm] to [yr]
    sig.yr <- 1e-2 * sig.cm * (HL$rho / bdot[i])
    
    # diffusion length in [yr] on equidistant time grid
    sigma <- stats::approx(t / t.res, sig.yr, seq_len(nt), rule = 2)$y

    return(sigma)

  })

  res <- data.frame(t.equi, sigma)

  colnames(res) <- c("Time", names)

  return(res)
  
}

#' Diffuse a record.
#' 
#' This function diffuses a time series or record with a given depth-dependent
#' diffusion length by convolution with a Gaussian kernel.
#'
#' This function expects a numeric vector with the depth-dependent diffusion
#' lengths of the same length as \code{rec}, or a single value to use a
#' constant diffusion length.
#'
#' The input diffusion length is internally scaled according to the resolution
#' of the record given by \code{res}. The convolution integral is then
#' solved by a simple summation over the kernel width set to ~ 10 times the
#' current diffusion length.
#'
#' For \code{debug = FALSE}: To avoid \code{NAs} at both ends of the diffused
#' version of \code{rec} resulting from the kernel extending beyond the record
#' ends, the kernel is clipped at the upper end to the range below the
#' surface. At the lower end, the record is extended by ~ 10 times the maximum
#' of \code{sigma} and filled with the mean average value of \code{rec}.
#'
#' @param rec Numeric vector containing the record that is to be diffused.
#' @param sigma Numeric vector of the diffusion lengths corresponding to the
#'   depths at which \code{rec} is tabulated, or of length one to diffuse
#'   \code{rec} with a constant diffusion length. Must be in units of the
#'   resolution of \code{rec} (typically [cm]).
#' @param res Resolution of \code{rec} in the same units as \code{sigma}.
#' @param debug if \code{TRUE} the values at top and bottom of the diffused
#'   record which are potentially affected by the finite record length are set
#'   to \code{NA}. Defaults to \code{FALSE}. See also Details.
#' @return Numeric vector containing the diffused version of \code{rec}.
#' @author Thomas Muench, modified by Thomas Laepple
DiffuseRecord <- function(rec, sigma, res = 1, debug = FALSE){

  if (missing(sigma)) {
    stop("No diffusion length passed as input.", call. = FALSE)
  }
  if (any(!is.finite(sigma))) {
    stop("Missing values passed as diffusion length.", call. = FALSE)
  }

  ns <- length(sigma)
  n <- length(rec)

  if (ns != 1 & ns != n) {
    stop("Diffusion length neither of length 1 nor matches length of record.",
         call. = FALSE)
  }

  # recycle sigma if needed
  if (ns == 1) {
    sigma <- rep(sigma, n)
  }

  # scale diffusion length according to resolution of record
  sigma <- sigma / res

  # pad end of record with mean of record to avoid NA's
  # at the end of diffused record
  if (!debug) rec <- c(rec, rep(mean(rec, na.rm = TRUE), 10 * max(sigma)))

  # vector to store diffused data
  rec.diffused <- rep(NA, n)

  # loop over record
  for (i in 1 : n) {

    # diffusion length for current depth
    sig <- sigma[i]

    if (sig == 0) {

      diff.value <- rec[i]

    } else {

      # set range of convolution integral (= 2*imax + 1) to ~ 10*sig
      imax <- ceiling(5 * sig)
      ran <- (i - imax) : (i + imax)

      # if part of range extends above surface, set diffused value to 'NA'
      # for 'debug = TRUE', else skip that part of range in the
      # convolution integral

      if (!all(ran > 0) & debug) {
        diff.value <- NA
      } else {
        ran <- ran[ran > 0]
        # relative range for convolution kernel
        rel.ran <- i - ran
        
        # convolution kernel
        kernel <- exp(-(rel.ran)^2 / (2 * sig^2))
        kernel <- kernel / sum(kernel)

        # diffuse data at current depth bin
        diff.value <- sum(rec[ran] * kernel)
      }
    }

    rec.diffused[i] <- diff.value
  }

  return(rec.diffused)

}
