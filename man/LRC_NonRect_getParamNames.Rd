% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRC_nonrectangular.R
\name{LRC_NonRect_getParamNames}
\alias{LRC_NonRect_getParamNames}
\title{\itemize{
\item \code{LRC_NonRect_getParamNames}: return the parameter names
used by this Light Response Curve Function
}}
\usage{
LRC_NonRect_getParamNames()
}
\value{
string vector of parameter names. Positions are important.
Adds sixth parameter, \code{logitconv} to the parameters of
\code{\link[=LRC_getParamNames]{LRC_getParamNames()}}
}
\description{
\itemize{
\item \code{LRC_NonRect_getParamNames}: return the parameter names
used by this Light Response Curve Function
}
}
\seealso{
\code{\link[=LRC_NonRect_predictGPP]{LRC_NonRect_predictGPP()}}
}
