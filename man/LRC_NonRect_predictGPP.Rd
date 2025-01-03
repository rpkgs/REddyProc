% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRC_nonrectangular.R
\name{LRC_NonRect_predictGPP}
\alias{LRC_NonRect_predictGPP}
\title{Nonrectangular Light Response function for GPP}
\usage{
LRC_NonRect_predictGPP(Rg, Amax, alpha, conv)
}
\arguments{
\item{Rg}{photosynthetic flux density \link{umol / m2 / s} or Global Radiation}

\item{Amax}{saturation (beta parameter) adjusted for effect of VPD}

\item{alpha}{slope at Rg = 0}

\item{conv}{convexity parameter (see details)}
}
\value{
vector of GPP
}
\description{
Nonrectangular Light Response function for GPP
}
\details{
This function generalizes the \code{\link[=RectangularLRCFitter_predictGPP]{RectangularLRCFitter_predictGPP()}}
by adding the convexity parameter \code{conv}.
For \code{conv -> 0 (logitconv -> -Inf)}: approaches the rectangular hyperbolic.
For \code{conv -> 1 (logitconv -> + Inf)}: approaches a step function.
Expected values of \code{conv} are about 0.7-0.9 (Moffat 2012).
}
\seealso{
\code{\link[=LRC_predictGPP]{LRC_predictGPP()}}
}
