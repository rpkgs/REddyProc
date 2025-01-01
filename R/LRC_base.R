#' R5 reference parent class for describing the NEP ~ PAR relationship
#'
#' @import methods
#' @export LightResponseCurveFitter
#' @exportClass LightResponseCurveFitter
LightResponseCurveFitter <- setRefClass("LightResponseCurveFitter")

# if this method is adjusted in subclass, then also adjust getPriorLocation,
# getPriorScale.

#' LightResponseCurveFitter_getParameterNames
#' 
#' @details 
#' - `k`     : VPD effect
#' - `beta`  : saturation of GPP at high radiation
#' - `alpha` : initial slope
#' - `RRef`  : basal respiration (units of provided NEE, usually mumol CO2 m-^-2 s^-2)
#' - `E0`    : temperature sensitivity estimated from night-time data (K)
#' 
#' @importFrom stats setNames
#' @export
LightResponseCurveFitter_getParameterNames <- function() {
  vars = c("k", "beta", "alpha", "RRef", "E0")
  setNames(vars, vars)
}

LightResponseCurveFitter$methods(
  getParameterNames = LightResponseCurveFitter_getParameterNames)

#' Optimize rectangular hyperbolic light response curve in one window
#' 
#' @details Optimization is performed for three initial parameter sets that
#' differ by beta0 (* 1.3, * 0.8). From those three, the optimization result is
#' selected that yielded the lowest misfit.
#' 
#' Starting values are: 
#' - `k` = 0, 
#' - `beta` = interpercentileRange(0.03, 0.97) of respiration, 
#' - `alpha` = 0.1, 
#' - `R_ref` from nightTime estimate.
#' - `E0` is fixed to the night-time estimate, but varies for estimating parameter uncertainty.
#' 
#' @param dsDay data.frame with columns NEE, Rg, Temp_C, VPD, and no NAs in NEE
#' @param E0 temperature sensitivity of respiration
#' @param sdE0 standard deviation of E_0.n
#' @param RRefNight basal respiration estimated from night time data
#' @param controlGLPart further default parameters (see \code{\link{partGLControl}})
#' @param lastGoodParameters numeric vector returned by last reasonable fit
#' 
#' @return A list with the following components:
#' - `thetaOpt`          : [numeric] vector of optimized parameters, including the fixed ones and E0
#' - `iOpt`              : [integer] vector of positions of parameters that are
#'   optimized, including E0, which has been optimized prior to this function
#' - `thetaInitialGuess` : [numeric] initial guess from data
#' - `covParms`          : [numeric] matrix of the covariance matrix of
#'   parameters, including E0
#' - convergence         : [integer] code specifying convergence problems
#'  + 0      : good convergence
#'  + 1-1000 : see \code{\link{optim}}
#'  + 1001   : too few bootstraps converged
#'  + 1002   : fitted parameters were outside reasonable bounds
#'  + 1003   : too few valid records in window
#'  + 1004   : near zero covariance in bootstrap indicating bad fit
#'  + 1005   : covariance from curvature of fit yielded negative variances
#'    indicating bad fit
#'  + 1006   : prediction of highest PAR in window was far from saturation
#'    indicating insufficient data to constrain LRC
#'  + 1010   : no temperature-respiration relationship found
#'  + 1011   : too few valid records in window (from different location:
#'    partGLFitLRCOneWindow)
#' 
#' @seealso [partGLFitLRCWindows()], [LightResponseCurveFitter_optimLRCBounds()]
#' @export
LightResponseCurveFitter_fitLRC <- function(
  dsDay, E0, sdE0, RRefNight, 
  controlGLPart = partGLControl(), lastGoodParameters = rep(NA_real_, 7L)) {

  # Three initial guess vectors are defined according to Lasslop et al., 2010
  parNames <- .self$getParameterNames() # hook method from derived classes
  nPar <- length(parNames)

  parPrior <- .self$getPriorLocation(dsDay$NEE, RRefNight = RRefNight, E0 = E0)
  thetaInitials <- .self$getParameterInitials(parPrior)
  
  resOpt3 <- apply(thetaInitials, 1, function(theta0) {
    resOpt <- .self$optimLRCBounds(theta0, parPrior,
      dsDay = dsDay, ctrl = controlGLPart, lastGoodParameters = lastGoodParameters)
  })
  iValid <- which(sapply(resOpt3, function(resOpt) { is.finite(resOpt$theta[1]) }))
  resOpt3Valid <- resOpt3[iValid]
  optSSE <- sapply(resOpt3Valid, "[[", "value")
  getNAResult <- function(convergenceCode) {
    list(
      thetaOpt = structure(rep(NA_real_, nPar), names = colnames(thetaInitials)),
      iOpt = integer(0), 
      thetaInitialGuess = thetaInitials[1, ],
      covParms = matrix(NA_real_, nPar, nPar, dimnames = list(
        colnames(thetaInitials), colnames(thetaInitials)
      )),
      convergence = convergenceCode
    )
  }
  if (sum(!is.na(optSSE)) == 0L) {
    return(getNAResult(resOpt3[[1]]$convergence))
  } else {
    resOpt <- resOpt3Valid[[iBest <- which.min(optSSE)]] # select lowest cost
    thetaOpt <- resOpt$theta
    if (controlGLPart$nBootUncertainty == 0L) {
      ## details<< If \code{controlGLPart$nBootUncertainty == 0L} then the
      ## covariance matrix of the
      ## parameters is estimated by the Hessian of the LRC curve at optimum.
      ## Then, the additional uncertainty and covariance with uncertainty E0
      ## is neglected.
      # seParmsHess <- seParmsHess0 <- sqrt(abs(diag(solve(resOpt$hessian))))
      covParmsLRC <- try(if ((resOpt$hessian[1L, 1L] < 1e-8)) {
        # case where k = 0 and not varying: cov(k, :) = cov(:, ) = 0
        covParmsLRC <- structure(diag(0, nrow = nrow(resOpt$hessian)),
          dimnames = dimnames(resOpt$hessian)
        )
        covParmsLRC[-1L, -1L] <- solve(resOpt$hessian[-1L, -1L])
        covParmsLRC
      } else {
        solve(resOpt$hessian)
      }, silent = TRUE)
      if (inherits(covParmsLRC, "try-error")) {
        return(getNAResult(1006L)) # count not invert the Hessian
      }
      covParms <- structure(diag(0, nrow = length(resOpt$theta)),
        dimnames = list(parNames, parNames)
      )
      covParms[5L, 5L] <- sdE0^2
      covParms[resOpt$iOpt, resOpt$iOpt] <- covParmsLRC
      if (any(diag(covParms) < 0)) {
        return(getNAResult(1005L))
      }
    } else {
      ## details<<
      ## If \code{controlGLPart.l$nBootUncertainty > 0L} then the
      ## covariance matrix of the
      ## parameters is estimated by a bootstrap of the data.
      ## In each draw, E0 is drawn from N ~ (E_0, sdE_0).
      # #seealso<< \code{\link{.bootStrapLRCFit}}
      resBoot <- .bootStrapLRCFit(resOpt$theta, resOpt$iOpt, dsDay, sdE0,
        parPrior, controlGLPart,
        lrcFitter = .self
      )
      iFiniteRows <- which(is.finite(resBoot[, 1L]))
      ## details<< If there are no estimates for more than 20% of the
      ## bootstrapped samples
      ## The an NA-result with convergence code 1001L is returned.
      if (length(iFiniteRows) < 0.8 * nrow(resBoot)) {
        return(getNAResult(1001L))
      }
      covParms <- cov(resBoot[iFiniteRows, ])
      if (covParms[2, 2] < 1e-8) getNAResult(1004L)
    }
  }
  # further parameter checking after parameter uncertainty has been computed
  # (before here, other parameter checking is done in optimLRCBounds)
  sdTheta <- thetaOpt
  sdTheta[] <- NA
  sdTheta[resOpt$iOpt] <- sqrt(diag(covParms)[resOpt$iOpt])
  if (!.self$isParameterInBounds(thetaOpt, sdTheta,
    RRefNight = RRefNight, ctrl = controlGLPart)) {
    return(getNAResult(1002L))
  }
  # debugging: tracing specific LRC fits:
  ## value<< a list, If none of the optimizations from different starting
  ## conditions converged, the parameters are NA.
  list(
    thetaOpt = thetaOpt, 
    iOpt = c(resOpt$iOpt, 5L), 
    thetaInitialGuess =  thetaInitials[1, ],
    covParms = covParms, 
    convergence = resOpt$convergence)
}
LightResponseCurveFitter$methods(fitLRC = LightResponseCurveFitter_fitLRC)


#' LightResponseCurveFitter_getPriorLocation
#' 
#' @param NEEDay numeric vector of daytime NEE
#' @param RRefNight numeric scalar of basal respiration estimated from night-time data
#' @param E0 numeric scalar of night-time estimate of temperature sensitivity
#' 
#' @return a numeric vector with prior estimates of the parameters
#' @export
LightResponseCurveFitter_getPriorLocation <- function(NEEDay, RRefNight, E0) {
  if (!is.finite(RRefNight)) stop("must provide finite RRefNight")
  zs = quantile(NEEDay, c(0.03, 0.97), na.rm = TRUE)
  c(k = 0.05, 
    beta = as.vector(abs(diff(zs))),
    alpha = 0.1, 
    RRef = as.vector(RRefNight), 
    E0 = as.vector(E0)
  ) # parameterPrior
}
LightResponseCurveFitter$methods(
  getPriorLocation = LightResponseCurveFitter_getPriorLocation)


#' LightResponseCurveFitter_getPriorScale
#' 
#' @details 
#' The beta parameter is quite well defined. Hence use a prior with
#' a standard deviation. The specific results are sometimes a bit sensitive
#' to the uncertainty of the beta prior. This uncertainty is set
#' corresponding to 20 times the median relative flux uncertainty. The prior
#' is weighted n times the observations in the cost. Hence, overall it is
#' using a weight of 1 / 20 of the weight of all observations.
#' 
#' However, its not well defined if PAR does not reach saturation. Need to
#' check before applying this prior
#' 
#' @param thetaPrior numeric vector of location of priors
#' @param medianRelFluxUncertainty numeric scalar: median across the
#' relative uncertainty of the flux values, i.e. sdNEE / NEE
#' @param nRec integer scalar: number of finite observations
#' @param ctrl list of further controls, with entry
#' `isLasslopPriorsApplied`
#' 
#' @return a numeric vector with prior estimates of the parameters
#' @export
LightResponseCurveFitter_getPriorScale <- function(
  thetaPrior, medianRelFluxUncertainty, nRec, ctrl) {
  sdParameterPrior <- if (ctrl$isLasslopPriorsApplied) {
    # twutz: changed to no prior for logitconv
    c(k = 50, beta = 600, alpha = 10, RRef = 80, E0 = NA)
  } else {
    sdBetaPrior <- 20 * medianRelFluxUncertainty * thetaPrior[2] / sqrt(nRec)
    c(k = NA, beta = as.vector(sdBetaPrior), alpha = NA, RRef = NA, E0 = NA)
  }
}
LightResponseCurveFitter$methods(getPriorScale = LightResponseCurveFitter_getPriorScale)

#' @export
LightResponseCurveFitter_getParameterInitials <- function(thetaPrior) {
  # only beta (second parameter)  is varied between different initial guesses
  ## value<< a numeric matrix (3, nPar) of initial values for fitting parameters
  theta0 <- matrix(rep(thetaPrior, each = 3), 3, length(thetaPrior),
    dimnames = list(NULL, names(thetaPrior))
  )
  theta0[2, 2] <- thetaPrior[2] * 1.3
  theta0[3, 2] <- thetaPrior[2] * 0.8
  theta0 # thetaInitials
}
LightResponseCurveFitter$methods(
  getParameterInitials = LightResponseCurveFitter_getParameterInitials)

#' Optimize parameters with refitting with some fixed parameters if outside bounds
#' 
#' @details 
#' If parameters alpha or k are outside bounds (Table A1 in Lasslop 2010),
#' refit with some parameters fixed to values from fit of previous window.
#' 
#' @param theta0 numeric vector of initial parameter estimate
#' @param parameterPrior numeric vector of prior estimate of model parameters
#' @param ... further parameters to \code{.optimLRC}
#' @param dsDay dataframe of NEE, sdNEE and predictors Rg, VPD and Temp
#' @param lastGoodParameters numeric vector of last successful fit
#' @param ctrl list of further controls, such as \code{isNeglectVPDEffect = TRUE}
#' 
#' @return list of result of [LightResponseCurveFitter_optimLRCOnAdjustedPrior()]
#' - `theta`: [numeric] vector of optimized parameters
#' - `iOpt`: [integer] vector of positions of parameters that are optimized
#' - `thetaInitialGuess`: [numeric] initial guess from data
#' see [LightResponseCurveFitter_fitLRC()]
#' 
#' @seealso [LightResponseCurveFitter_fitLRC()]
#' @export
LightResponseCurveFitter_optimLRCBounds <- function(
  theta0, parameterPrior, ..., dsDay, lastGoodParameters, ctrl) {
  # twutz 161014: default alpha
  if (!is.finite(lastGoodParameters[3L])) lastGoodParameters[3L] <- 0.22
  isNeglectVPDEffect <- isTRUE(ctrl$isNeglectVPDEffect)
  VPD0 <- 10 # VPD0 fixed to 10 hPa, Lasslop 2010

  isUsingFixedVPD <- isNeglectVPDEffect || (sum(dsDay$VPD >= VPD0, na.rm = TRUE) == 0)
  isUsingFixedAlpha <- FALSE
  
  getIOpt <- .self$getOptimizedParameterPositions
  
  theta0Adj <- theta0 # initial estimate with some parameters adjusted to bounds
  if (isNeglectVPDEffect) theta0Adj[1] <- 0
  resOpt <- resOpt0 <- .self$optimLRCOnAdjustedPrior(theta0Adj,
    iOpt = getIOpt(isUsingFixedVPD, isUsingFixedAlpha),
    parameterPrior = parameterPrior, ctrl, dsDay = dsDay, ...)
  
  fun_optim <- function() {
    .self$optimLRCOnAdjustedPrior(theta0Adj,
      iOpt = getIOpt(isUsingFixedVPD, isUsingFixedAlpha),
      parameterPrior = parameterPrior, ctrl, dsDay = dsDay, ...
    )
  }

  # dsDay <- list(ctrl, ...)$dsDay
  if (is.na(resOpt$theta[1L]) || (resOpt$theta[1L] < 0)) { # k
    isUsingFixedVPD <- TRUE
    theta0Adj[1L] <- 0
    resOpt <- fun_optim()
    # check alpha, in case refit with fixed alpha of last window
    if ((is.na(resOpt$theta[3L]) || (resOpt$theta[3L] > 0.22)) &&
      is.finite(lastGoodParameters[3L])) {
      isUsingFixedAlpha <- TRUE
      theta0Adj[3L] <- lastGoodParameters[3L]
      resOpt <- fun_optim()
    }
  } else {
    # check alpha, if gt 0.22 estimate parameters with fixed alpha of last window
    # if not last window exists, let alpha > 0.22
    if ((is.na(resOpt$theta[3L]) || (resOpt$theta[3L] > 0.22)) &&
      is.finite(lastGoodParameters[3L])) {
      isUsingFixedAlpha <- TRUE
      theta0Adj[3L] <- lastGoodParameters[3L]
      resOpt <- fun_optim()
      # check k, if less than zero estimate parameters without VPD effect
      # and with fixed alpha of last window
      if (is.na(resOpt$theta[1L]) || (resOpt$theta[1L] < 0)) {
        isUsingFixedVPD <- TRUE
        theta0Adj[1L] <- 0
        resOpt <- fun_optim()
      }
    }
  }
  ## details<<
  ## No parameters are reported if alpha<0 or RRef < 0 or beta0 < 0
  ## or beta0 > 250
  # positions in theta0: "k"     "beta0" "alpha"  "RRef"    "E0"
  if (resOpt$convergence != 0) {
    resOpt$theta <- NA
  }
  if (!is.na(resOpt$theta[1L]) && ((resOpt$theta[3L] < 0) ||
    (resOpt$theta[4L] < 0) || (resOpt$theta[2L] < 0) ||
    (resOpt$theta[2L] >= 250))
  ) {
    # TODO estimate RRef from daytime data?
    # LloydT_E0fix
    # stop("case with alpha or beta < 0")
    resOpt$theta[] <- NA
    resOpt$convergence <- 1002
  }
  ## details<<
  ## Not parameters are reported if the data did not contain records that
  ## are near light saturation.
  ## This is checked by comparing the prediction at highest PAR with the
  ## beta parameter
  if (is.finite(resOpt$theta[1]) & !isTRUE(ctrl$isUsingLasslopQualityConstraints) &
    is.finite(ctrl$minPropSaturation)
  ) {
    dsDay <- list(...)$dsDay
    iMaxRg <- which.max(dsDay$Rg)
    dsDayMax <- dsDay[iMaxRg, , drop = FALSE]
    # compute prediction at maximum observed PAR and and PAR = 2000 and
    # compare how close its to saturation
    predMaxGPP <- .self$predictLRC(
      theta = resOpt$theta, 
      Rg = c(dsDayMax$Rg, 2000),
      VPD = 0 # dsDayMax$VPD
      , Temp = NA # dsDayMax$Temp
      # , VPD0 = 10 			# TODO: think of providing VPD0 and TRef to function
      # , TRef = 15
    )$GPP
    if (predMaxGPP[1] < ctrl$minPropSaturation * predMaxGPP[2]) {
      # plot(-NEE ~ Rg, dsDay)
      resOpt$theta[] <- NA
      resOpt$convergence <- 1006
    }
  }
  # Further checks are done, after parameter uncertainty has been determined,
  # by call from fitLRC to isParameterInBounds
  resOpt
}
LightResponseCurveFitter$methods(
  optimLRCBounds = LightResponseCurveFitter_optimLRCBounds)

#' Get the positions of the parameters to optimize for given Fixed
#' 
#' If subclasses extend the parameter vector, they need to override this method.
#' 
#' @param isUsingFixedVPD boolean scalar: if TRUE, VPD effect set to zero and is not optimized
#' @param isUsingFixedAlpha boolean scalar: if TRUE, initial slope is fixed and is not optimized
#' 
#' @return integer vector of positions in parameter vector
#' @export
LightResponseCurveFitter_getOptimizedParameterPositions <- function(isUsingFixedVPD, isUsingFixedAlpha) {
  if (!isUsingFixedVPD & !isUsingFixedAlpha) {
    c(1:4)
  } else if (isUsingFixedVPD & !isUsingFixedAlpha) {
    2:4
  } else if (!isUsingFixedVPD & isUsingFixedAlpha) {
    c(1L, 2L, 4L)
  } else if (isUsingFixedVPD & isUsingFixedAlpha) {
    c(2L, 4L)
  } # iOpt
}
LightResponseCurveFitter$methods(
  getOptimizedParameterPositions = LightResponseCurveFitter_getOptimizedParameterPositions)

#' optimLRCOnAdjustedPrior
#' 
#' Lower bound flux uncertainty and adjust prior uncertainty before calling optimLRC
#' 
#' @details Only those records are used for optimization where both NEE and sdNEE are finite.
#' In larger settings, already filtered at \code{partGLFitLRCOneWindow}
#' 
#' Optimization of LRC parameters takes into account the uncertainty of the flux
#' values. In order to avoid very strong leverage, values with a very low
#' uncertainty (< a lower quantile) are assigned the lower quantile is assigned.
#' This procedure downweighs records with a high uncertainty, but does not apply
#' a large leverage for records with a very low uncertainty. Avoid this
#' correction by setting \code{ctrl$isBoundLowerNEEUncertainty = FALSE}
#' 
#' @param theta numeric vector of starting values
#' @param iOpt integer vector: positions of subset of parameters that are optimized
#' @param dsDay dataframe of NEE, sdNEE and predictors Rg, VPD and Temp
#' @param parameterPrior numeric vector of prior parameter estimates (corresponding to theta)
#' @param ctrl list of further controls
#' @param ... further arguments to
#' \code{\link{LightResponseCurveFitter_optimLRC}} (passed to
#' \code{\link{LightResponseCurveFitter_computeCost}})
#' 
#' @return list of result of [LightResponseCurveFitter_optimLRC()] amended with list
#' `theta`, `iOpt` and `convergence`
#' @export
LightResponseCurveFitter_optimLRCOnAdjustedPrior <- function(
  theta, iOpt, dsDay, parameterPrior, ctrl, ...) {
  
  if (!all(is.finite(theta))) stop("need to provide finite starting values.")
  dsDayFinite <- dsDay[is.finite(dsDay$NEE) & is.finite(dsDay$sdNEE), ]
  if (nrow(dsDayFinite) < ctrl$minNRecInDayWindow) {
    stop(
      "inspect too few records, should be already filtered ",
      "in partGLFitLRCOneWindow")
    return(list(
      theta = {
        theta[] <- NA
        theta
      }, iOpt = integer(0), convergence = 1003L))
  }

  minUnc <- quantile(dsDayFinite$sdNEE, 0.3)
  if (minUnc == 0) {
    stop("Too many zeros in uncertainty of NEE.",
      " This cannot be handled in daytime partitioning.")
  }
  Fc_unc <- if (isTRUE(ctrl$isBoundLowerNEEUncertainty)) {
    # twutz: avoid excessive weights by small uncertainties (of 1 / unc^2)
    pmax(dsDayFinite$sdNEE, minUnc)
  } else {
    dsDayFinite$sdNEE
  }
  # plot(Fc_unc ~ dsDayFinite$sdNEE)
  medianRelFluxUncertainty <- abs(median(Fc_unc / dsDayFinite$NEE))
  ## details<<
  ## The uncertainty of the prior, that maybe derived from fluxes)  is allowed to
  ## adapt to the uncertainty of the fluxes.
  ## This is done in \code{link{LightResponseCurveFitter_getPriorScale}}
  sdParameterPrior <- .self$getPriorScale(parameterPrior, medianRelFluxUncertainty,
    nrow(dsDayFinite),
    ctrl = ctrl
  )
  sdParameterPrior[-iOpt] <- NA
  isUsingHessian <- (ctrl$nBootUncertainty == 0L)
  .self$optimLRC(theta,
    iOpt = iOpt,
    flux = -dsDayFinite$NEE,
    sdFlux = Fc_unc,
    sdParameterPrior = sdParameterPrior,
    parameterPrior = parameterPrior,
    Rg = dsDayFinite$Rg,
    VPD = dsDayFinite$VPD,
    Temp = dsDayFinite$Temp,
    isUsingHessian = isUsingHessian,
    ctrl = ctrl
  )
}

LightResponseCurveFitter$methods(
  optimLRCOnAdjustedPrior = LightResponseCurveFitter_optimLRCOnAdjustedPrior)

#' isParameterInBounds
#' 
#' Check if estimated parameter vector is within reasonable bounds.
#' Check the Beta bounds that depend on uncertainty: outside if (beta > 100 and
#' sdBeta >= beta)
#' 
#' @param theta numeric vector of parameters
#' @param sdTheta numeric vector of standard deviation of parameters
#' @param RRefNight numeric scalar: night-time based estimate of basal respiration
#' @param ctrl list of further controls
#' 
#' @return logical scalar: TRUE if parameters are within bounds, FALSE otherwise
#' @export
LightResponseCurveFitter_isParameterInBounds <- function(theta, sdTheta, RRefNight, ctrl) {
  if (!is.finite(theta[2])) return(FALSE)
  if (isTRUE(as.vector((theta[2] > 100) && (sdTheta[2] >= theta[2])))) return(FALSE)
  return(TRUE)
}
LightResponseCurveFitter$methods(
  isParameterInBounds = LightResponseCurveFitter_isParameterInBounds)

#' optim LRC
#' 
#' @param theta numeric vector of parameters
#' @param iOpt integer vector of positions of parameters that are optimized
#' @param sdParameterPrior numeric vector of standard deviation of parameterPrior
#' @param ... other arguments to \code{\link{LightResponseCurveFitter_computeCost}}
#' @param ctrl list of further controls
#' @param isUsingHessian scalar boolean: set to TRUE to compute Hessian at optimum
#' 
#' @return list of result of \code{\link{optim}} amended with list
#' @export
LightResponseCurveFitter_optimLRC <- function(
  theta, iOpt, sdParameterPrior, ..., ctrl, isUsingHessian) {
  thetaOrig <- theta
  resOptim <- optim(thetaOrig[iOpt], .self$computeCost,
    theta = thetaOrig, iOpt = iOpt,
    sdParameterPrior = sdParameterPrior, ...,
    control = list(reltol = ctrl$LRCFitConvergenceTolerance),
    method = "BFGS", hessian = isUsingHessian
  )
  
  thetaOpt <- theta
  thetaOpt[iOpt] <- resOptim$par
  ans <- list(theta = thetaOpt, iOpt = iOpt)
  c(resOptim, ans)
}
LightResponseCurveFitter$methods(optimLRC = LightResponseCurveFitter_optimLRC)


#' Computing residual sum of squares for predictions vs. data of NEE
#' 
#' @param thetaOpt numeric vector of optimized parameters
#' @param theta numeric vector of parameters
#' @param iOpt integer vector of positions of parameters that are optimized
#' @param flux numeric: NEP (-NEE) or GPP time series [umolCO2 / m2 / s],
#' should not contain NA
#' @param sdFlux numeric: standard deviation of Flux [umolCO2 / m2 / s],
#' should not contain NA
#' @param parameterPrior numeric vector along theta: prior estimate of
#' parameter (range of values)
#' @param sdParameterPrior standard deviation of parameterPrior
#' @param weightMisfitPar2000 weight of misfit of difference between
#' saturation and prediction at PAR = 2000
#' @param ... other arguments to
#' \code{\link{LightResponseCurveFitter_predictLRC}}, such as VPD0, fixVPD
#' @return numeric: residual sum of squares
#' 
#' @export
LightResponseCurveFitter_computeCost <- function(
  thetaOpt, theta, iOpt, flux, sdFlux, parameterPrior, sdParameterPrior, ...) 
{
  theta[iOpt] <- thetaOpt
  resPred <- .self$predictLRC(theta, ...)
  NEP_mod <- resPred$NEP
  # if (is.na(mean(NEP_mod)) == TRUE) {
  #  recover()
  # }
  misFitPrior <- (((theta - parameterPrior)) / (sdParameterPrior))^2
  misFitObs <- sum(((NEP_mod - flux) / sdFlux)^2)
  RSS <- misFitObs + sum(misFitPrior, na.rm = TRUE)
  # if (!is.finite(RSS) ) recover()	# debugging the fit
  RSS
}
LightResponseCurveFitter$methods(
  computeCost = LightResponseCurveFitter_computeCost)

#' Light Response Function
#'
#' Predicts the Net Ecosystem Production (NEP), Ecosystem Respiration (Reco) and
#' Gross Primary Production (GPP) based on the Light Response Function.
#'
#' The VPD effect is included according to Lasslop et al., 2010.
#'
#' @param theta numeric vector of parameters. If theta is a matrix, a different
#' row of parameters is used for different entries of other inputs.
#' @param Rg [numeric] -> photosynthetic flux density [mumol / m2 / s] or
#' Global Radiation.
#' @param VPD [numeric] -> Vapor Pressure Deficit [hPa]
#' @param Temp [numeric] -> Temperature [degC]
#' @param VPD0 [hPa] -> Parameters VPD0 fixed to 10 hPa according to Lasslop et al 2010
#' @param fixVPD if TRUE the VPD effect is not considered and VPD is not part of the computation
#' @param TRef numeric scalar of Temperature (degree Celsius) for reference respiration RRef
#' 
#' @return
#' - NEP: Net ecosystem production (-NEE), vector of length(Rg)
#' - Reco: Ecosystem respiration
#' - GPP: Gross primary production
#'
#' @export
LightResponseCurveFitter_predictLRC <- function(
    theta, Rg, VPD, Temp, VPD0 = 10, fixVPD = (k == 0), TRef = 15) {
  if (is.matrix(theta)) {
    k <- theta[, 1]
    beta <- theta[, 2]
    alpha <- theta[, 3]
    RRef <- theta[, 4]
    E0 <- theta[, 5]
  } else {
    k <- theta[1]
    beta <- theta[2]
    alpha <- theta[3]
    RRef <- theta[4]
    E0 <- theta[5]
  }
  if (length(fixVPD) != length(VPD)) {
    if (length(fixVPD) == 1L) {
      fixVPD <- rep(fixVPD, length(VPD))
    } else {
      stop("Length of vector argument fixVPD must correspond to rows in theta.")
    }
  }
  Amax <- ifelse(fixVPD, beta,
    # ifelse(is.finite(VPD) & (VPD > VPD0), beta * exp(-k * (VPD-VPD0)), beta)
    # deprecated: better filter for k: twutz: 170927: introduced pmin(1, ...)
    # after looking at pvWave code, can happen if k is negative
    ifelse((VPD > VPD0), beta * exp(-k * (VPD - VPD0)), beta)
  )
  Reco <- RRef * exp(E0 * (1 / ((273.15 + TRef) - 227.13) - 1 /
    (Temp + 273.15 - 227.13)))
  GPP <- .self$predictGPP(Rg, Amax = Amax, alpha = alpha)
  NEP <- GPP - Reco
  ## a data.frame of length of Rg of computed
  ans <- list(NEP = NEP, Reco = Reco, GPP = GPP)
}
LightResponseCurveFitter$methods(predictLRC = LightResponseCurveFitter_predictLRC)

#' Light Response function for GPP
#'
#' - Rectangular: [RectangularLRCFitter_predictGPP()]
#' - Nonrectangular: [NonrectangularLRCFitter_predictGPP()]
#' - Logistic Sigmoid: [LogisticSigmoidLRCFitter_predictGPP()]
#'
#' @param Rg ppfd [numeric] -> photosynthetic flux density [mumol / m2 / s] or Global Radiation
#' @param ... further parameters to the LRC
#' @return numeric vector of length(Rg) of GPP
#'
#' @seealso [partitionNEEGL()]
#' @export
LightResponseCurveFitter_predictGPP <- function(Rg, ...) {
  stop("Abstract method. Need to define in derived LRC class.")
}
LightResponseCurveFitter$methods(predictGPP = LightResponseCurveFitter_predictGPP)


#' Gradient of the Light Response Function
#' @param theta numeric vector of parameters. If theta is a matrix, a different
#' row of parameters is used for different entries of other inputs.
#' @param Rg ppfd [numeric] -> photosynthetic flux density [mumol / m2 / s] or Global Radiation
#' @param VPD VPD [numeric] -> Vapor Pressure Deficit [hPa]
#' @param Temp Temp [numeric] -> Temperature [degC]
#' @param VPD0 VPD0 [hPa] -> Parameters VPD0 fixed to 10 hPa according to Lasslop et al 2010
#' @param fixVPD if TRUE the VPD effect is not considered and VPD is not part of the computation
#' @param TRef numeric scalar of Temperature (degree Celsius) for reference respiration RRef
#' 
#' @return list with gradient matrices. For each record (length(Rg)), c("k", "beta", "alpha", "RRef")
#' @export
LightResponseCurveFitter_computeLRCGradient <- function(
  theta, Rg, VPD, Temp, VPD0 = 10, fixVPD = (k == 0), TRef = 15) {
  if (is.matrix(theta)) {
    k <- theta[, 1]
    beta <- theta[, 2]
    alpha <- theta[, 3]
    RRef <- theta[, 4]
    E0 <- theta[, 5]
  } else {
    k <- theta[1]
    beta <- theta[2]
    alpha <- theta[3]
    RRef <- theta[4]
    E0 <- theta[5]
  }
  if (length(fixVPD) != length(VPD)) {
    if (length(fixVPD) == 1L) {
      fixVPD <- rep(fixVPD, length(VPD))
    } else {
      stop("Length of vector argument fixVPD must correspond to rows in theta.")
    }
  }
  Amax <- ifelse(fixVPD, beta,
    # ifelse(is.finite(VPD) & (VPD > VPD0), beta * exp(-k * (VPD-VPD0)), beta)
    ifelse((VPD > VPD0), beta * exp(-k * (VPD - VPD0)), beta)
  )
  # ex <- expression(beta * exp(-k * (VPD - VPD0)) ); deriv(ex, c("beta", "k"))
  dAmax_dkVPD <- ifelse(fixVPD, 0,
    ifelse(VPD > VPD0, beta * -(VPD - VPD0) * exp(-k * (VPD - VPD0)), 0)
  )
  dAmax_dbeta0 <- ifelse(fixVPD, 0,
    ifelse(VPD > VPD0, exp(-k * (VPD - VPD0)), 1)
  )
  # Reco<- RRef * exp(E0 * (1 / ((273.15 + 10)-227.13)-1 / (Temp + 273.15-227.13)))
  # ex <- expression(RRef * exp(E0 * (1 / ((273.15 + TRef)-227.13)-1 /
  #  (Temp + 273.15-227.13))) ); deriv(ex, c("RRef", "E0"))
  .expr7 <- 1 / (273.15 + TRef - 227.13) - 1 / (Temp + 273.15 - 227.13)
  .expr9 <- exp(E0 * .expr7)
  gradReco <- matrix(0,
    ncol = 2L, nrow = length(.expr9), dimnames =
      list(NULL, c("RRef", "E0"))
  )
  gradReco[, "RRef"] <- dReco_dRRef <- .expr9
  gradReco[, "E0"] <- dReco_dE0 <- RRef * (.expr9 * .expr7)
  
  gradGPP <- array(0, c(nrow(gradReco), 3L), list(NULL, c("k", "beta", "alpha")))
  dGPP_dAMax <- .self$computeGPPGradient(Rg, Amax, alpha)
  gradGPP[, "beta"] <- dGPP_dAMax[, 1] * dAmax_dbeta0
  gradGPP[, "k"] <- dGPP_dAMax[, 1] * dAmax_dkVPD
  gradGPP[, "alpha"] <- dGPP_dAMax[, 2]
  # NEP <- GPP - Reco
  gradNEP <- cbind(gradGPP, -gradReco)
  ## list with gradient matrices. For each record
  ## (length(Rg)), c("k", "beta", "alpha", "RRef")
  list(NEP = gradNEP, Reco = gradReco, GPP = gradGPP)
}

LightResponseCurveFitter$methods(
  computeLRCGradient = LightResponseCurveFitter_computeLRCGradient)
