# random, percentile, density and quantile function of the logit-normal distribution
# and estimation of parameters from percentiles by Sum of Squares Newton optimization
#

### Transforming (0,1) to normal scale (-Inf Inf)
## details<<
## function \eqn{ logit(p)= log \left( \frac{p}{1-p} \right) = log(p) - log(1-p) }
## seealso<< \code{\link{invlogit}}
logit <- function(p, ...) {
  qlogis(p, ...)
}

### Transforming (-Inf,Inf) to original scale (0,1)
## details<<
## function \eqn{f(z) = \frac{e^{z}}{e^{z} + 1} \! = \frac{1}{1 + e^{-z}} \!}
## seealso<< \code{\link{logit}}
invlogit <- function(q, ...) {
  plogis(q, ...)
}

#' replace NA by value of previous record
#'
#' @param x numeric vector to fill NAs
#' @param firstValue value to be used for NA at the beginning of x
#'
#' @return numeric vector with NAs replaced
fillNAForward <- function(x, firstValue = median(x, na.rm = TRUE)) {
  iMissing <- which(!is.finite(x))
  if (length(iMissing) && (iMissing[1] == 1L)) {
    x[1L] <- firstValue
    iMissing <- iMissing[-1]
  }
  if (length(iMissing)) {
    for (i in iMissing) {
      # set to value from previous window
      x[i] <- x[i - 1L]
    }
  }
  return(x)
}

#' replace missing standard deviation of a measure x by a percentage of x
#'
#' @details If either perc or inSdX is NA then only the other criterion is applied.
#' If both are NA then all missings are set to NA.
#'
#' @param sdX numeric vector: with missing to be repalce
#' @param x numeric vector of length(sdX): value form which percentage is computed
#' @param perc numeric scalar: sdX = perc * x
#' @param minSdX numeric scalar: minimum of sdX to be applied for low x
#'
#' @return sdX with non-finite values replaced.
replaceMissingSdByPercentage <- function(sdX, x, perc = 0.2, minSdX = 0.7) {
  ## \code{sdX[iToFill] <- pmax(minSdX, abs(x[iToFill] * perc), na.rm = TRUE)}
  iToFill <- !is.finite(sdX)
  sdX[iToFill] <- pmax(minSdX, abs(x[iToFill] * perc), na.rm = TRUE)
  sdX
}


listk <- function(...) {
  # get variable names from input expressions
  cols <- as.list(substitute(list(...)))[-1]
  vars <- names(cols)
  Id_NoName <- if (is.null(vars)) seq_along(cols) else which(vars == "")

  if (length(Id_NoName) > 0) {
    vars[Id_NoName] <- sapply(cols[Id_NoName], deparse)
  }
  setNames(list(...), vars)
}
