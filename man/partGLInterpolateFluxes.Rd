% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PartitioningLasslop10.R
\name{partGLInterpolateFluxes}
\alias{partGLInterpolateFluxes}
\title{Interpolate Reco and GPP from two neighboring parameter sets}
\usage{
partGLInterpolateFluxes(
  Rg,
  VPD,
  Temp,
  resParms,
  controlGLPart = partGLControl(),
  lrcFitter,
  isVerbose = TRUE
)
}
\arguments{
\item{Rg}{numeric vector of photosynthetic flux density \link{umol / m2 / s} or Global Radiation}

\item{VPD}{numeric vector of Vapor Pressure Deficit \link{hPa}}

\item{Temp}{numeric vector of Temperature \link{degC}}

\item{resParms}{data frame with results of \code{partGLFitLRCWindows} of
fitting the light-response-curve for several windows}

\item{controlGLPart}{further default parameters, see \code{\link{partGLControl}}}

\item{lrcFitter}{R5 class instance responsible for fitting the LRC curve}

\item{isVerbose}{set to FALSE to suppress messages}
}
\description{
Interpolate Reco and GPP from two neighboring parameter sets
}
\details{
\code{resLRC$iFirstRecInCentralDay} must denote the row for which the
LRC parameters are representative, here, the first record of the center day
}
\seealso{
\code{\link[=partitionNEEGL]{partitionNEEGL()}}, \code{\link[=partGLFitLRCWindows]{partGLFitLRCWindows()}}
}
