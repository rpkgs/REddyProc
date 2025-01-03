% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRC_base_optim.R
\name{LRC_predictLRC}
\alias{LRC_predictLRC}
\title{Light Response Function}
\usage{
LRC_predictLRC(theta, Rg, VPD, Temp, VPD0 = 10, fixVPD = (k == 0), TRef = 15)
}
\arguments{
\item{theta}{numeric vector of parameters. If theta is a matrix, a different
row of parameters is used for different entries of other inputs.}

\item{Rg}{photosynthetic flux density \link{mumol / m2 / s} or
Global Radiation.}

\item{VPD}{Vapor Pressure Deficit \link{hPa}}

\item{Temp}{Temperature \link{degC}}

\item{VPD0}{\link{hPa}, Parameters VPD0 fixed to 10 hPa according to Lasslop et al 2010}

\item{fixVPD}{if TRUE the VPD effect is not considered and VPD is not part of the computation}

\item{TRef}{numeric scalar of Temperature (degree Celsius) for reference respiration RRef}
}
\value{
\itemize{
\item NEP: Net ecosystem production (-NEE), vector of length(Rg)
\item Reco: Ecosystem respiration
\item GPP: Gross primary production
}
}
\description{
Predicts the Net Ecosystem Production (NEP), Ecosystem Respiration (Reco) and
Gross Primary Production (GPP) based on the Light Response Function.
}
\details{
The VPD effect is included according to Lasslop et al., 2010.
}
