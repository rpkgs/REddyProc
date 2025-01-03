% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PartitioningLasslop10.R
\name{partGLFitLRCOneWindow}
\alias{partGLFitLRCOneWindow}
\title{Estimate Rectangular Hyperbolic Light Response Curve parameters}
\usage{
partGLFitLRCOneWindow(
  ds,
  winInfo,
  prevRes,
  E0Win,
  controlGLPart = partGLControl(),
  lrcFitter = RectangularLRCFitter()
)
}
\arguments{
\item{ds}{data.frame with numeric columns NEE, sdNEE, Temp (degC), VPD, Rg,
and logical columns isNight and isDay}

\item{winInfo}{one-row data.frame with window information, including iWindow}

\item{prevRes}{component prevRes from previous result, here with item prevE0}

\item{E0Win}{data.frame with columns E0, sdE0, RRef from nighttime,
one row for each window}

\item{controlGLPart}{list of further default parameters}

\item{lrcFitter}{R5 class instance responsible for fitting the light response curve}
}
\value{
a list with components resOpt, summary, isValid
}
\description{
Estimate Rectangular Hyperbolic Light Response Curve parameters
}
\details{
Estimation as in Lasslop et al., 2010 for successive periods, i.e.
windows. Fitting is done on a subset of the data where isDay, NEE, sdNEE,
Temp, Rg, and VPD are all non-NA and isDay is TRUE
}
\seealso{
\code{\link[=partitionNEEGL]{partitionNEEGL()}}
}
