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
