#' optim LRC
#'
#' @param theta numeric vector of parameters
#' @param iOpt integer vector of positions of parameters that are optimized
#' @param sdParameterPrior numeric vector of standard deviation of parameterPrior
#' @param ... other arguments to [LRC_computeCost()]
#' @param ctrl list of further controls
#' @param isUsingHessian scalar boolean: set to TRUE to compute Hessian at optimum
#'
#' @return list of result of [optim()] amended with list
#' @export
LRC_optimLRC <- function(
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
LightResponseCurveFitter$methods(optimLRC = LRC_optimLRC)


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
#' [LRC_predictLRC()], such as VPD0, fixVPD
#' @return numeric: residual sum of squares
#'
#' @export
LRC_computeCost <- function(
    thetaOpt, theta, iOpt, flux, sdFlux, parameterPrior, sdParameterPrior, ...) {
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
  computeCost = LRC_computeCost
)

#' Light Response Function
#'
#' Predicts the Net Ecosystem Production (NEP), Ecosystem Respiration (Reco) and
#' Gross Primary Production (GPP) based on the Light Response Function.
#'
#' The VPD effect is included according to Lasslop et al., 2010.
#'
#' @param theta numeric vector of parameters. If theta is a matrix, a different
#' row of parameters is used for different entries of other inputs.
#' @param Rg -> photosynthetic flux density [mumol / m2 / s] or
#' Global Radiation.
#' @param VPD -> Vapor Pressure Deficit [hPa]
#' @param Temp -> Temperature [degC]
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
LRC_predictLRC <- function(
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
  Reco <- cal_Reco(Temp, E0, RRef, TRef)
  GPP <- .self$predictGPP(Rg, Amax = Amax, alpha = alpha)
  NEP <- GPP - Reco
  ## a data.frame of length of Rg of computed
  ans <- list(NEP = NEP, Reco = Reco, GPP = GPP)
}
LightResponseCurveFitter$methods(predictLRC = LRC_predictLRC)

cal_Reco <- function(Temp, E0, RRef, TRef = 15) {
  RRef * exp(E0 * (1 / ((273.15 + TRef) - 227.13) - 1 / (Temp + 273.15 - 227.13)))
}

#' Light Response function for GPP
#'
#' - Rectangular: [RectangularLRCFitter_predictGPP()]
#' - Nonrectangular: [NonrectangularLRCFitter_predictGPP()]
#' - Logistic Sigmoid: [LogisticSigmoidLRCFitter_predictGPP()]
#'
#' @param Rg ppfd -> photosynthetic flux density [mumol / m2 / s] or Global Radiation
#' @param ... further parameters to the LRC
#' @return numeric vector of length(Rg) of GPP
#'
#' @seealso [partitionNEEGL()]
#' @export
LRC_predictGPP <- function(Rg, ...) {
  stop("Abstract method. Need to define in derived LRC class.")
}
LightResponseCurveFitter$methods(predictGPP = LRC_predictGPP)


#' Gradient of the Light Response Function
#' @param theta numeric vector of parameters. If theta is a matrix, a different
#' row of parameters is used for different entries of other inputs.
#' @param Rg ppfd -> photosynthetic flux density [mumol / m2 / s] or Global Radiation
#' @param VPD VPD -> Vapor Pressure Deficit [hPa]
#' @param Temp Temp -> Temperature [degC]
#' @param VPD0 VPD0 [hPa] -> Parameters VPD0 fixed to 10 hPa according to Lasslop et al 2010
#' @param fixVPD if TRUE the VPD effect is not considered and VPD is not part of the computation
#' @param TRef numeric scalar of Temperature (degree Celsius) for reference respiration RRef
#'
#' @return list with gradient matrices. For each record (length(Rg)), c("k", "beta", "alpha", "RRef")
#' @export
LRC_computeLRCGradient <- function(
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
  computeLRCGradient = LRC_computeLRCGradient
)
