#' R5 reference parent class for describing the NEP ~ PAR relationship
#'
#' @import methods
#' @export LightResponseCurveFitter
#' @exportClass LightResponseCurveFitter
LightResponseCurveFitter <- setRefClass("LightResponseCurveFitter")

# if this method is adjusted in subclass, then also need to adjust getPriorLocation,
# getPriorScale.

#' LRC_getParameterNames
#' 
#' @details 
#' - `k`     : VPD effect
#' - `beta`  : saturation of GPP at high radiation
#' - `alpha` : initial slope
#' - `RRef`  : basal respiration (units of provided NEE, usually mumol CO2 m-2 s-2)
#' - `E0`    : temperature sensitivity estimated from night-time data (K)
#' 
#' @importFrom stats setNames
#' @export
LRC_getParameterNames <- function() {
  vars = c("k", "beta", "alpha", "RRef", "E0")
  setNames(vars, vars)
}
LightResponseCurveFitter$methods(
  getParameterNames = LRC_getParameterNames)

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
#' @param controlGLPart further default parameters (see [partGLControl()])
#' @param lastGoodParameters numeric vector returned by last reasonable fit
#' 
#' @return A list with the following components:
#' - `thetaOpt`          : vector of optimized parameters, including the fixed ones and E0
#' - `iOpt`              : vector of positions of parameters that are optimized,
#'   including E0, which has been optimized prior to this function
#' - `thetaInitialGuess` : initial guess from data
#' - `covParms`          : matrix of the covariance matrix of parameters,
#'   including E0
#' - convergence         : code specifying convergence problems
#'  + 0      : good convergence
#'  + 1-1000 : see [optim()]
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
#' @seealso [partGLFitLRCWindows()], [LRC_optimLRCBounds()]
#' @export
LRC_fitLRC <- function(
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
  iValid <- which(sapply(resOpt3, \(r) is.finite(r$theta[1]) ))
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
        parPrior, controlGLPart, lrcFitter = .self)
      iFiniteRows <- which(is.finite(resBoot[, 1L]))
      ## details<< If there are no estimates for more than 20% of the
      ## bootstrapped samples
      ## The an NA-result with convergence code 1001L is returned.
      if (length(iFiniteRows) < 0.8 * nrow(resBoot)) return(getNAResult(1001L))
      
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
LightResponseCurveFitter$methods(fitLRC = LRC_fitLRC)


#' LRC_getPriorScale
#' 
#' @details 
#' The beta parameter is quite well defined. Hence use a prior with
#' a standard deviation. The specific results are sometimes a bit sensitive
#' to the uncertainty of the beta prior. 
#' 
#' This uncertainty is set corresponding to 20 times the median relative flux
#' uncertainty. The prior is weighted n times the observations in the cost.
#' Hence, overall it is using a weight of 1 / 20 of the weight of all
#' observations.
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
LRC_getPriorScale <- function(
  thetaPrior, medianRelFluxUncertainty, nRec, ctrl) {
  sdParameterPrior <- if (ctrl$isLasslopPriorsApplied) {
    # twutz: changed to no prior for logitconv
    c(k = 50, beta = 600, alpha = 10, RRef = 80, E0 = NA)
  } else {
    sdBetaPrior <- 20 * medianRelFluxUncertainty * thetaPrior[2] / sqrt(nRec)
    c(k = NA, beta = as.vector(sdBetaPrior), alpha = NA, RRef = NA, E0 = NA)
  }
}
LightResponseCurveFitter$methods(getPriorScale = LRC_getPriorScale)

#' LRC_getPriorLocation
#' 
#' @param NEEDay numeric vector of daytime NEE
#' @param RRefNight numeric scalar of basal respiration estimated from night-time data
#' @param E0 numeric scalar of night-time estimate of temperature sensitivity
#'
#' @return a numeric vector with prior estimates of the parameters
#' @export
LRC_getPriorLocation <- function(NEEDay, RRefNight, E0) {
  if (!is.finite(RRefNight)) stop("must provide finite RRefNight")
  zs <- quantile(NEEDay, c(0.03, 0.97), na.rm = TRUE)
  c(
    k = 0.05,
    beta = as.vector(abs(diff(zs))),
    alpha = 0.1,
    RRef = as.vector(RRefNight),
    E0 = as.vector(E0)
  ) # parameterPrior
}
LightResponseCurveFitter$methods(getPriorLocation = LRC_getPriorLocation)

# ' @return A numeric matrix (3, nPar) of initial values for fitting parameters
#' @export
LRC_getParameterInitials <- function(thetaPrior) {
  theta0 <- matrix(rep(thetaPrior, each = 3), 3, length(thetaPrior),
    dimnames = list(NULL, names(thetaPrior))
  )
  theta0[2, 2] <- thetaPrior[2] * 1.3 # jitter `beta`
  theta0[3, 2] <- thetaPrior[2] * 0.8
  theta0 # thetaInitials
}
LightResponseCurveFitter$methods(
  getParameterInitials = LRC_getParameterInitials)

#' Optimize parameters with refitting with some fixed parameters if outside bounds
#' 
#' @details 
#' If parameters alpha or k are outside bounds (Table A1 in Lasslop 2010),
#' refit with some parameters fixed to values from fit of previous window.
#' 
#' @param theta0 numeric vector of initial parameter estimate
#' @param parameterPrior numeric vector of prior estimate of model parameters
#' @param ... further parameters to `.optimLRC`
#' @param dsDay dataframe of NEE, sdNEE and predictors Rg, VPD and Temp
#' @param lastGoodParameters numeric vector of last successful fit
#' @param ctrl list of further controls, such as `isNeglectVPDEffect = TRUE`
#' 
#' @return list of result of [LRC_optimLRCOnAdjustedPrior()]
#' - `theta`: vector of optimized parameters
#' - `iOpt`: vector of positions of parameters that are optimized
#' - `thetaInitialGuess`: initial guess from data
#' see [LRC_fitLRC()]
#' 
#' @seealso [LRC_fitLRC()]
#' @export
LRC_optimLRCBounds <- function(
  theta0, parameterPrior, ..., dsDay, lastGoodParameters, ctrl) {
  # twutz 161014: default alpha
  if (!is.finite(lastGoodParameters[3L])) lastGoodParameters[3L] <- 0.22 # default alpha
  
  isNeglectVPDEffect <- isTRUE(ctrl$isNeglectVPDEffect)
  VPD0 <- 10 # VPD0 fixed to 10 hPa, Lasslop 2010

  isUsingFixedVPD <- isNeglectVPDEffect || (sum(dsDay$VPD >= VPD0, na.rm = TRUE) == 0)
  isUsingFixedAlpha <- FALSE
  
  getIOpt <- .self$getOptimizedParameterPositions
  
  theta0Adj <- theta0 # initial estimate with some parameters adjusted to bounds
  if (isNeglectVPDEffect) theta0Adj[1] <- 0
  resOpt <- resOpt0 <- .self$optimLRCOnAdjustedPrior(theta0Adj,
    iOpt = getIOpt(isUsingFixedVPD, isUsingFixedAlpha),
    parameterPrior = parameterPrior, ctrl, dsDay = dsDay, ...
  )
  
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
      , Temp = NA # dsDayMax$Temp #! BUG here, Temp is not provided
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
  optimLRCBounds = LRC_optimLRCBounds)

#' Get the positions of the parameters to optimize for given Fixed
#' 
#' If subclasses extend the parameter vector, they need to override this method.
#' 
#' @param isUsingFixedVPD boolean scalar: if TRUE, VPD effect set to zero and is not optimized
#' @param isUsingFixedAlpha boolean scalar: if TRUE, initial slope is fixed and is not optimized
#' 
#' @return integer vector of positions in parameter vector
#' @export
LRC_getOptimizedParameterPositions <- function(isUsingFixedVPD, isUsingFixedAlpha) {
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
  getOptimizedParameterPositions = LRC_getOptimizedParameterPositions)

#' optimLRCOnAdjustedPrior
#' 
#' Lower bound flux uncertainty and adjust prior uncertainty before calling optimLRC
#' 
#' @details Only those records are used for optimization where both NEE and sdNEE are finite.
#' In larger settings, already filtered at `partGLFitLRCOneWindow`
#' 
#' Optimization of LRC parameters takes into account the uncertainty of the flux
#' values. In order to avoid very strong leverage, values with a very low
#' uncertainty (< a lower quantile) are assigned the lower quantile is assigned.
#' This procedure downweighs records with a high uncertainty, but does not apply
#' a large leverage for records with a very low uncertainty. Avoid this
#' correction by setting `ctrl$isBoundLowerNEEUncertainty = FALSE`
#' 
#' @param theta numeric vector of starting values
#' @param iOpt integer vector: positions of subset of parameters that are optimized
#' @param dsDay dataframe of NEE, sdNEE and predictors Rg, VPD and Temp
#' @param parameterPrior numeric vector of prior parameter estimates (corresponding to theta)
#' @param ctrl list of further controls
#' @param ... further arguments to
#' [LRC_optimLRC()] (passed to
#' [LRC_computeCost()])
#' 
#' @return list of result of [LRC_optimLRC()] amended with list
#' `theta`, `iOpt` and `convergence`
#' @export
LRC_optimLRCOnAdjustedPrior <- function(
  theta, iOpt, dsDay, parameterPrior, ctrl, ...) {
  
  if (!all(is.finite(theta))) stop("need to provide finite starting values.")
  dsDayFinite <- subset(dsDay, is.finite(NEE) & is.finite(sdNEE))

  npar = length(theta)
  if (nrow(dsDayFinite) < ctrl$minNRecInDayWindow) {
    stop("inspect too few records, should be already filtered ",
      "in partGLFitLRCOneWindow")
    return(list(
      theta = rep(NA, npar), iOpt = integer(0), convergence = 1003L))
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
  
  medianRelFluxUncertainty <- abs(median(Fc_unc / dsDayFinite$NEE))
  ## details<<
  ## The uncertainty of the prior, that maybe derived from fluxes)  is allowed to
  ## adapt to the uncertainty of the fluxes.
  ## This is done in [LRC_getPriorScale()]
  sdParameterPrior <- .self$getPriorScale(parameterPrior, medianRelFluxUncertainty,
    nrow(dsDayFinite), ctrl = ctrl)
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
  optimLRCOnAdjustedPrior = LRC_optimLRCOnAdjustedPrior)

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
LRC_isParameterInBounds <- function(theta, sdTheta, RRefNight, ctrl) {
  if (!is.finite(theta[2])) return(FALSE)
  if (isTRUE(as.vector((theta[2] > 100) && (sdTheta[2] >= theta[2])))) return(FALSE)
  return(TRUE)
}
LightResponseCurveFitter$methods(
  isParameterInBounds = LRC_isParameterInBounds)
