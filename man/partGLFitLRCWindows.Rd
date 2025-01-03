% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PartitioningLasslop10.R
\name{partGLFitLRCWindows}
\alias{partGLFitLRCWindows}
\title{Estimate successive Rectangular Hyperbolic Light Response Curve parameters}
\usage{
partGLFitLRCWindows(
  ds,
  winSizeDays = 4L,
  strideInDays = 2L,
  nRecInDay = 48L,
  dsTempSens,
  isVerbose = TRUE,
  controlGLPart = partGLControl(),
  lrcFitter = RectangularLRCFitter()
)
}
\arguments{
\item{ds}{data.frame with numeric columns NEE, sdNEE, Temp (degC), VPD, Rg,
and logical columns isNight and isDay}

\item{winSizeDays}{Window size in days for daytime fits}

\item{strideInDays}{step in days for shifting the windows}

\item{nRecInDay}{number of records within one day (for half-hourly data its 48)}

\item{dsTempSens}{data.frame that reports for each window temperature sensitivity
parameters E0 and RRef}

\item{isVerbose}{set to FALSE to suppress messages}

\item{controlGLPart}{list of further default parameters}

\item{lrcFitter}{R5 class instance responsible for fitting the light response curve}
}
\value{
a tibble with a row for each window providing information on window size,
parameters estimated from light response curve (LRC), and their standard deviation.
Moreover, and estimated from night-time data after smoothing and forward-filling:
the uncertainty of temperature sensitivity \code{E0_night_sd}, And the respiration
at reference temperature \code{RRef_night}. Finally it contains a list column
\code{resOpt} with the optimization results.
}
\description{
Estimate successive Rectangular Hyperbolic Light Response Curve parameters
}
\seealso{
\code{\link[=partGLFitLRCOneWindow]{partGLFitLRCOneWindow()}}
}
