#' Logistic sigmoid Light response curve
#' 
#' @import methods
#' @export LogisticSigmoidLRCFitter
#' @exportClass LogisticSigmoidLRCFitter
LogisticSigmoidLRCFitter <- setRefClass("LogisticSigmoidLRCFitter",
  contains = "LightResponseCurveFitter")

#' Logistic Sigmoid Light Response function for GPP
#'
#' @details `GPP <- Amax * tanh(alpha * Rg / Amax)`
#'
#' @param Rg ppfd -> photosynthetic flux density [mumol / m2 / s] or Global Radiation
#' @param Amax saturation (beta parameter) adjusted for effect of VPD
#' @param alpha slope at Rg = 0
#'
#' @seealso [LRC_predictGPP()]
#' @return
#' - `GPP`: vector of GPP
#' - `GPPGradient`: [numeric matrix] of gradients of predicted GPP to Amax and alpha
#'
#' @rdname LogisticSigmoidLRCFitter
#' @export
LogisticSigmoidLRCFitter_predictGPP <- function(Rg, Amax, alpha) {
  Amax * tanh(alpha * Rg / Amax) # GPP
}

LogisticSigmoidLRCFitter$methods(predictGPP = LogisticSigmoidLRCFitter_predictGPP)

#' @rdname LogisticSigmoidLRCFitter
#' @export
LogisticSigmoidLRCFitter_computeGPPGradient <- function(Rg, Amax, alpha) {
  # ex <- expression(  Amax * tanh(alpha * Rg / Amax) ); deriv(ex, c("Amax", "alpha"))
  .expr1 <- alpha * Rg
  .expr2 <- .expr1 / Amax
  .expr3 <- tanh(.expr2)
  .expr8 <- cosh(.expr2)^2
  .value <- Amax * .expr3
  .grad <- array(0, c(length(.value), 2L), list(NULL, c("Amax", "alpha")))
  .grad[, 1L] <- .expr3 - Amax * (.expr1 / Amax^2 / .expr8)
  .grad[, 2L] <- Amax * (Rg / Amax / .expr8)
  .grad
}

LogisticSigmoidLRCFitter$methods(
  computeGPPGradient = LogisticSigmoidLRCFitter_computeGPPGradient)
