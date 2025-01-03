% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PartitioningLasslop10Nighttime.R
\name{partGL_FitNight_1win_RRef}
\alias{partGL_FitNight_1win_RRef}
\title{Estimate Reference temperature from nighttime and given tempsens E0}
\usage{
partGL_FitNight_1win_RRef(
  dss,
  winInfo,
  prevRes = list(),
  E0Win,
  controlGLPart = partGLControl(),
  TRef = 15
)
}
\arguments{
\item{dss}{data.frame with numeric columns NEE, isNight, Temp, Rg}

\item{winInfo}{one-row data.frame with window information, including iWindow}

\item{prevRes}{component prevRes from previous result, here not used.}

\item{E0Win}{data.frame with columns E0 and sdE0, RRefFit, and TRefFit
with one row for each window}

\item{controlGLPart}{list of further default parameters}

\item{TRef}{numeric scalar of Temperature (degree Celsius) for reference
respiration RRef}
}
\value{
named scalar, RRef
}
\description{
Estimation of respiration at reference temperature (RRef) and temperature
sensitivity (E0) for one window.
}
\details{
If there are too few records
(n < \code{controlGLPart$minNRecInDayWindow}) then return NA.
}
\seealso{
\code{\link[=partitionNEEGL]{partitionNEEGL()}}, \code{\link[=partGL_FitNight_1win_E0_kernel]{partGL_FitNight_1win_E0_kernel()}}
}
