#' Nonrectangular Light response curve
#'
#' @import methods
#' @export NonrectangularLRCFitter
#' @exportClass NonrectangularLRCFitter
NonrectangularLRCFitter <- setRefClass("NonrectangularLRCFitter",
  contains = "LightResponseCurveFitter"
)

#' - `LRC_NonRect_getParamNames`: return the parameter names
#'   used by this Light Response Curve Function
#'
#' @return string vector of parameter names. Positions are important.
#' Adds sixth parameter, `logitconv` to the parameters of
#' [LRC_getParamNames()]
#'
#' @seealso [LRC_NonRect_predictGPP()]
#' @export
LRC_NonRect_getParamNames <- function() {
  ans <- callSuper()
  c(ans, logitconf = "logitconv") ## << logit-transformed convexity parameter.
  ## The value at original scale is obtained by
  ## \code{conv = 1 / (1 + exp(-logitconv))}
}

NonrectangularLRCFitter$methods(
  getParamNames = LRC_NonRect_getParamNames)

# LRC_getPriorLocation
LRC_NonRect_getPriorLocation <- function(NEEDay, RRefNight, E0) {
  ans <- callSuper(NEEDay = NEEDay, RRefNight = RRefNight, E0 = E0)
  c(ans, logitconv = logit(0.75))
}
NonrectangularLRCFitter$methods(getPriorLocation = LRC_NonRect_getPriorLocation)

#' Get the prior distribution of parameters
#' @keywords internal
#' @inheritParams LRC_getPriorScale
LRC_NonRect_getPriorScale <- function(thetaPrior, medianRelFluxUncertainty, nRec, ctrl) {
  ans <- callSuper(thetaPrior, medianRelFluxUncertainty, nRec, ctrl)
  c(ans, logitconv = NA)
}
NonrectangularLRCFitter$methods(getPriorScale = LRC_NonRect_getPriorScale)

#' Get the positions of the parameters to optimize for given Fixed
#'
#' @param isUsingFixedVPD boolean: if TRUE, VPD effect set to zero and is not optimized
#' @param isUsingFixedAlpha boolean: if TRUE, initial slope is fixed and is not optimized
#'
#' @return integer vector of positions of the parameters to optimize
LRC_NonRect_getOptimizedParamPos <- function(
    isUsingFixedVPD, isUsingFixedAlpha) {
  iOpt <- callSuper(isUsingFixedVPD = isUsingFixedVPD, isUsingFixedAlpha = isUsingFixedAlpha)
  c(iOpt, 6) # add the convexity parameter
}
NonrectangularLRCFitter$methods(
  getOptimizedParamPos =
    LRC_NonRect_getOptimizedParamPos
)

#' Nonrectangular Hyperbolic Light Response function: (Gilmanov et al., 2003)
#'
#' @param theta numeric vector of parameters
#' @param Rg photosynthetic flux density [umol / m2 / s] or Global Radiation
#' @param VPD Vapor Pressure Deficit [hPa]
#' @param Temp Temperature [degC]
#' @param VPD0 [hPa] -> Parameters VPD0 fixed to 10 hPa according to
#' Lasslop et al 2010
#' @param fixVPD boolean scalar or vector of nrow theta: fixVPD if TRUE the VPD
#' effect is not considered and VPD is not part of the computation
#' @param TRef numeric scalar of Temperature (degree Celsius) for reference respiration RRef
#' 
#' @export
LRC_NonRect_predictLRC <- function(
    theta, Rg, VPD, Temp, VPD0 = 10, fixVPD = (k == 0), TRef = 15) {
  if (is.matrix(theta)) {
    k <- theta[, 1]
    beta <- theta[, 2]
    alpha <- theta[, 3]
    RRef <- theta[, 4]
    E0 <- theta[, 5]
    logitconv <- theta[, 6]
  } else {
    k <- theta[1]
    beta <- theta[2]
    alpha <- theta[3]
    RRef <- theta[4]
    E0 <- theta[5]
    logitconv <- theta[6]
  }

  conv <- invlogit(logitconv)
  if (length(fixVPD) != length(VPD)) {
    if (length(fixVPD) == 1L) {
      fixVPD <- rep(fixVPD, length(VPD))
    } else {
      stop("Length of vector argument fixVPD must correspond to rows in theta.")
    }
  }
  Amax <- ifelse(fixVPD, beta,
    ifelse((VPD > VPD0), beta * exp(-k * (VPD - VPD0)), beta))
  Reco <- RRef * exp(E0 * (1 / ((273.15 + TRef) - 227.13)
    - 1 / (Temp + 273.15 - 227.13)))
  GPP <- .self$predictGPP(Rg, Amax = Amax, alpha = alpha, conv = conv)
  NEP <- GPP - Reco
  list(NEP = NEP, Reco = Reco, GPP = GPP)
}
NonrectangularLRCFitter$methods(predictLRC = LRC_NonRect_predictLRC)

#' Nonrectangular Light Response function for GPP
#' 
#' @details 
#' This function generalizes the [RectangularLRCFitter_predictGPP()]
#' by adding the convexity parameter `conv`.
#' For `conv -> 0 (logitconv -> -Inf)`: approaches the rectangular hyperbolic.
#' For `conv -> 1 (logitconv -> + Inf)`: approaches a step function.
#' Expected values of `conv` are about 0.7-0.9 (Moffat 2012).
#' 
#' @param Rg photosynthetic flux density [umol / m2 / s] or Global Radiation
#' @param Amax saturation (beta parameter) adjusted for effect of VPD
#' @param alpha slope at Rg = 0
#' @param conv convexity parameter (see details)
#' 
#' @seealso [LRC_predictGPP()]
#' @return vector of GPP
#' @export
LRC_NonRect_predictGPP <- function(Rg, Amax, alpha, conv) {
  zRoot <- ((alpha * Rg + Amax)^2) - (4 * alpha * Rg * conv * Amax)
  zRoot[which(zRoot < 0)] <- 0
  GPP <- (1 / (2 * conv)) * (alpha * Rg + Amax - sqrt(zRoot))
}
NonrectangularLRCFitter$methods(predictGPP = LRC_NonRect_predictGPP)


#' Gradient of [LRC_NonRect_predictLRC()]
#' 
#' @details This function computes the gradient of the Nonrectangular Light Response
#' Curve function with respect to the parameters.
#' 
#' Differs from base by extracting `conv` parameter from `theta` and adding
#' gradient to `logitconv` (3rd parameter from computeGPPGradient).
#' 
#' @param theta -> parameter vector (
#' theta[1] = kVPD (k),
#' theta[2] = beta0 (beta),
#' theta[3] = alpha,
#' theta[4] = RRef (rb), # theta[4] = E0,
#' theta[5] = logitconv)
#' @inheritParams LRC_predictLRC
#' 
#' @export
LRC_NonRect_computeLRCGradient <- function(
    theta, Rg, VPD, Temp, VPD0 = 10, fixVPD = (k == 0), TRef = 15) {
  if (is.matrix(theta)) {
    k <- theta[, 1]
    beta <- theta[, 2]
    alpha <- theta[, 3]
    RRef <- theta[, 4]
    E0 <- theta[, 5]
    logitconv <- theta[, 6]
  } else {
    k <- theta[1]
    beta <- theta[2]
    alpha <- theta[3]
    RRef <- theta[4]
    E0 <- theta[5]
    logitconv <- theta[6]
  }
  if (!is.finite(logitconv[1])) stop("need to provide finite logitconv in theta")
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
  # ex <- expression(beta * exp(-k * (VPD-VPD0)) ); deriv(ex, c("beta", "k"))
  dAmax_dkVPD <- ifelse(fixVPD, 0,
    ifelse(VPD > VPD0, beta * -(VPD - VPD0) * exp(-k * (VPD - VPD0)), 0))
  dAmax_dbeta0 <- ifelse(fixVPD, 0,
    ifelse(VPD > VPD0, exp(-k * (VPD - VPD0)), 1))
  # Reco<- RRef * exp(E0 * (1 / ((273.15 + 10)-227.13)-1 / (Temp + 273.15-227.13)))
  # ex <- expression(RRef * exp(E0 * (1 / ((273.15 + TRef)-227.13)-1 / (Temp + 273.15-227.13))) ); deriv(ex, c("RRef", "E0"))
  # to prevent numeric instabilities, use do not let temperature go below 20degC
  # .expr7 <- 1 / (273.15 + TRef - 227.13) - 1 / (Temp + 273.15 - 227.13)
  .expr7 <- 1 / (273.15 + TRef - 227.13) - 1 / (pmax(-20, Temp) + 273.15 - 227.13)
  .expr9 <- exp(E0 * .expr7)
  gradReco <- matrix(0,
    ncol = 2L, nrow = length(.expr9), dimnames =
      list(NULL, c("RRef", "E0"))
  )
  gradReco[, "RRef"] <- dReco_dRRef <- .expr9
  gradReco[, "E0"] <- dReco_dE0 <- RRef * (.expr9 * .expr7)
  
  gradGPP <- array(0, c(nrow(gradReco), 4L), list(NULL, c("k", "beta", "alpha", "logitconv")))
  dGPP_dAMax <- .self$computeGPPGradient(Rg, Amax, alpha, logitconv)
  gradGPP[, "beta"] <- dGPP_dAMax[, 1L] * dAmax_dbeta0
  gradGPP[, "k"] <- dGPP_dAMax[, 1L] * dAmax_dkVPD
  gradGPP[, "alpha"] <- dGPP_dAMax[, 2L]
  gradGPP[, "logitconv"] <- dGPP_dAMax[, 3L]
  # NEP <- GPP - Reco
  gradNEP <- cbind(gradGPP, -gradReco)
  ## value<<
  ## list with gradient matrices. For each record
  ## (length(Rg)), c("k", "beta", "alpha", "RRef")
  list(NEP = gradNEP, Reco = gradReco, GPP = gradGPP)
}

NonrectangularLRCFitter$methods(
  computeLRCGradient = LRC_NonRect_computeLRCGradient)

.tmp.f <- function() {
  iNonFinite <- which(!is.finite(gradReco[, "E0"]))
}

#' Gradient of [LRC_NonRect_predictGPP()]
#' 
#' @inheritParams LRC_NonRect_predictGPP
#' @param logitconv -> logit-transformed convexity parameter
LRC_NonRect_computeGPPGradient <- function(Rg, Amax, alpha, logitconv) {
  zRoot <- ((alpha * Rg + Amax)^2) - (4 * alpha * Rg * invlogit(logitconv) * Amax)
  iNegRoot <- which(zRoot < 0)
  # GPP<- (1 / (2 * (1 / (1 + exp(-logitconv))) )) * (alpha * Rg + Amax-sqrt(zRoot))
  # GPP<-               (1 / (2 * (1 / (1 + exp(-logitconv))) )) * (alpha * Rg + Amax-sqrt(((alpha * Rg + Amax)^2)-(4 * alpha * Rg * (1 / (1 + exp(-logitconv))) * Amax)))
  # ex <- expression(  (1 / (2 * (1 / (1 + exp(-logitconv))) )) * (alpha * Rg + Amax-sqrt(((alpha * Rg + Amax)^2)-(4 * alpha * Rg * (1 / (1 + exp(-logitconv))) * Amax))) ); deriv(ex, c("Amax", "alpha", "logitconv"))
  .expr2 <- exp(-logitconv)
  .expr3 <- 1 + .expr2
  .expr4 <- 1 / .expr3
  .expr5 <- 2 * .expr4
  .expr6 <- 1 / .expr5
  .expr8 <- alpha * Rg + Amax
  .expr11 <- 4 * alpha * Rg
  .expr12 <- .expr11 * .expr4
  .expr14 <- .expr8^2 - .expr12 * Amax
  .expr16 <- .expr8 - sqrt(.expr14)
  .expr20 <- .expr14^-0.5
  .expr36 <- .expr2 / .expr3^2
  .value <- .expr6 * .expr16
  # plot(.value ~ GPP)
  .grad <- array(0, c(length(.value), 3L), list(NULL, c("Amax", "alpha", "logitconv")))
  .grad[, "Amax"] <- .expr6 * (1 - 0.5 * ((2 * .expr8 - .expr12) * .expr20))
  .grad[, "alpha"] <- .expr6 * (Rg - 0.5 * ((2 * (Rg * .expr8) -
    4 * Rg * .expr4 * Amax) * .expr20))
  .grad[, "logitconv"] <- .expr6 * (0.5 * (.expr11 * .expr36 *
    Amax * .expr20)) - 2 * .expr36 / .expr5^2 * .expr16
  if (length(iNegRoot)) {
    # GPP<- (1 / (2 * (1 / (1 + exp(-logitconv))) )) * (alpha * Rg + Amax-0)
    # ex <- expression(  (1 / (2 * (1 / (1 + exp(-logitconv))) )) * (alpha * Rg + Amax-0) ); deriv(ex, c("Amax", "alpha", "logitconv"))
    # .expr2 <- exp(-logitconv)
    # .expr3 <- 1 + .expr2
    .expr5 <- 2 * (1 / .expr3)
    .expr6 <- 1 / .expr5
    .expr9 <- alpha * Rg + Amax - 0
    .value <- .expr6 * .expr9
    .gradNegRoot <- array(0, c(length(.value), 3L), list(NULL, c(
      "Amax",
      "alpha", "logitconv"
    )))
    .gradNegRoot[, "Amax"] <- .expr6
    .gradNegRoot[, "alpha"] <- .expr6 * Rg
    .gradNegRoot[, "logitconv"] <- -(2 * (.expr2 / .expr3^2) / .expr5^2 * .expr9)
    .grad[iNegRoot, ] <- .gradNegRoot[iNegRoot, ]
  }
  .grad
}

NonrectangularLRCFitter$methods(
  computeGPPGradient = LRC_NonRect_computeGPPGradient)
