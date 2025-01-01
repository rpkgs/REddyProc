#' R5 reference class for the Rectangular Light response curve
#'
#' @import methods
#' @export RectangularLRCFitter
#' @exportClass RectangularLRCFitter
RectangularLRCFitter <- setRefClass('RectangularLRCFitter'
                                    , contains = 'LightResponseCurveFitter')

#' RectangularLRCFitter_predictGPP
#' 
#' @details `GPP <- (Amax * alpha * Rg) / (alpha * Rg + Amax)`
#' 
#' @param Rg ppfd -> photosynthetic flux density [mumol / m2 / s] or Global Radiation
#' @param Amax saturation (beta parameter) adjusted for effect of VPD
#' @param alpha slope at Rg = 0
#' 
#' @return numeric vector of GPP
#' @seealso [LRC_predictGPP()]
#' @export
RectangularLRCFitter_predictGPP <- function(Rg, Amax, alpha) {
	(Amax * alpha * Rg) / (alpha * Rg + Amax) # GPP
}
RectangularLRCFitter$methods(predictGPP = RectangularLRCFitter_predictGPP)


#' RectangularLRCFitter_computeGPPGradient
#' 
#' @inheritParams RectangularLRCFitter_predictGPP
#' 
#' @return numeric matrix (length(Rg), 2) of gradients of predicted GPP
#' to Amax and alpha
RectangularLRCFitter_computeGPPGradient <- function(Rg, Amax, alpha) {
  # ex <- expression( (Amax * alpha * Rg) / (alpha * Rg + Amax) );
  # deriv(ex, c("Amax", "alpha"))
  .expr2 <- Amax * alpha * Rg
  .expr3 <- alpha * Rg
  .expr4 <- .expr3 + Amax
  .expr7 <- .expr4^2
  .value <- .expr2 / .expr4
  .grad <- array(0, c(length(.value), 2L), list(NULL, c("Amax", "alpha")))
  .grad[, 1L] <- .expr3 / .expr4 - .expr2 / .expr7
  .grad[, 2L] <- Amax * Rg / .expr4 - .expr2 * Rg / .expr7
  .grad
}
RectangularLRCFitter$methods(
  computeGPPGradient = RectangularLRCFitter_computeGPPGradient)

#' R5 reference class C-version for the Rectangular Light response curve
#'
#' @import methods
#' @export RectangularLRCFitterCVersion
#' @exportClass RectangularLRCFitterCVersion
RectangularLRCFitterCVersion <- setRefClass('RectangularLRCFitterCVersion',
                                            contains = 'RectangularLRCFitter'
	### overiding computeCost of \code{\link{RectangularLRCFitter}}
	### to a C-version (\code{\link{RectangularLRCFitter_C_computeCost}}).
)

#' Computing residual sum of squares for predictions vs. data of NEE
#' 
#' @param thetaOpt parameter vector with components of theta0 that are optimized
#' @param theta parameter vector with positions as in argument of
#' [LRC_getParamNames()]
#' @param iOpt position in theta that are optimized
#' @param flux numeric: NEP (-NEE) or GPP time series [umolCO2 / m2 / s],
#' should not contain NA
#' @param sdFlux numeric: standard deviation of Flux [umolCO2 / m2 / s],
#' should not contain NA
#' @param parameterPrior numeric vector along theta: prior estimate of
#' parameter (range of values)
#' @param sdParameterPrior standard deviation of parameterPrior
#' @param ... other arguments to
#' [LRC_predictLRC], such as VPD0, fixVPD
#' @param VPD0 VPD0 [hPa] -> Parameters VPD0 fixed to 10 hPa according
#' to Lasslop et al 2010
#' @param fixVPD boolean scalar or vector of nrow theta:fixVPD
#' if TRUE the VPD effect is not considered and VPD is not
#' part of the computation
#' 
#' @export
RectangularLRCFitterCVersion_computeCost <- function(
    thetaOpt, theta, iOpt, flux, sdFlux, parameterPrior, sdParameterPrior, ...,
    VPD0 = 10, fixVPD = (k == 0)) {
  theta[iOpt] <- thetaOpt
  k <- theta[1] # here k and fixVPD is only a scalar
  RHLightResponseCostC(
    theta, flux, sdFlux, parameterPrior, sdParameterPrior, ...,
    VPD0 = VPD0, fixVPD = fixVPD)
}

RectangularLRCFitterCVersion$methods(
  computeCost = RectangularLRCFitterCVersion_computeCost)
