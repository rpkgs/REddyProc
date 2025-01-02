#' compute logical vector of each rows in ds is its a valid night record
#' @details For robustness, data is trimmed to conditions at temperature > 1 degC
#' but only timmed if there are more at least 12 records left
#' @param ds data.frame with columns isNight, NEE, Temp (degC)
#' @return logical vector of length nrow(ds)
isValidNightRecord <- function(ds) {
  isValid <- !is.na(ds$isNight) & ds$isNight & !is.na(ds$NEE) & is.finite(ds$Temp)
  isFreezing <- ds$Temp[isValid] <= -1
  if (sum(!isFreezing) >= 12L) isValid[isValid][isFreezing] <- FALSE
  return(isValid)
}

#' Estimate temperature sensitivity E0 within one window
#' 
#' @details 
#' Estimation of respiration at reference temperature (RRef) and temperature
#' sensitivity (E0) for one window of night-time data.
#' 
#' The reference temperature is taken as the median of the temperatures of
#' valid records, unless a fixed reference temperature (in degree Celsius) is
#' specified by argument `controlGLPart$fixedTRefAtNightTime`.
#' 
#' @param dss data.frame with numeric columns NEE, sdNEE, Temp (degC) , VPD, Rg,
#' and logical columns isNight and isDay
#' @param winInfo one-row data.frame with window information, including iWindow
#' @param prevRes component prevRes from previous result, here with item prevE0
#' @param isVerbose set to FALSE to suppress messages
#' @param nRecInDay number of records within one day (for half-hourly data its 48)
#' @param controlGLPart list of further default parameters
#' 
#' @seealso [partitionNEEGL()], [partGLEstimateTempSensInBoundsE0Only()]
partGLFitNightTempSensOneWindow <- function(
  dss, winInfo, prevRes, isVerbose = TRUE, nRecInDay = 48L, 
  controlGLPart = partGLControl()) 
{
  isValid <- isValidNightRecord(dss)
  # check that there are enough night and enough day-values for fitting,
  # else continue with next window
  if (sum(isValid) < controlGLPart$minNRecInDayWindow) 
    return(data.frame(E0 = NA_real_, sdE0 = NA_real_, TRefFit = NA_real_, RRefFit = NA_real_))
  
  dssNight <- dss[isValid, ]

  TRefFitCelsius <- if (length(controlGLPart$fixedTRefAtNightTime) &&
    is.finite(controlGLPart$fixedTRefAtNightTime)
  ) {
    controlGLPart$fixedTRefAtNightTime
  } else {
    median(dssNight$Temp, na.rm = TRUE)
  }
  resNightFit <- partGLEstimateTempSensInBoundsE0Only(
    dssNight$NEE,
    fConvertCtoK(dssNight$Temp),
    prevE0 = prevRes$E0, TRefFit = fConvertCtoK(TRefFitCelsius))
  as.data.frame(resNightFit) # E0, sdE0, TRefFit, RRefFit
}


#' Estimate bounded temperature sensitivity E0 and RRef of ecosystem respiration
#' 
#' @details Basal respiration is reported for temperature of 15 degree Celsius.
#' However during the fit a reference temperature of the median of the dataset
#' is used. This is done to avoid strong correlations between estimated
#' parameters `E0` and `RRef`, that occure if reference temperature is outside
#' the center of the data.
#' 
#' @param REco night time NEE
#' @param temperatureKelvin temperature in Kelvin
#' @param prevE0 numeric scalar: the previous guess of Temperature Sensitivity
#' @param TRefFit numeric scalar reference temperature for which RRef is estimated.
#' Set it to the center of the data (the default) for best fit (lowest uncertainty)
#' 
#' @seealso [partGLFitLRCWindows()]
#' 
#' @return list with entries
#' - `E0`      : estimated temperature sensitivty E0 bounded to [50, 400]
#' - `sdE0`    : standard deviation of E0
#' - `TRefFit` : reference temperature used in the E0 fit
#' - `RRefFit` : respiration at TRefFit
#' @export 
partGLEstimateTempSensInBoundsE0Only <- function(
  REco, temperatureKelvin, prevE0 = NA, 
  TRefFit = median(temperatureKelvin, na.rm = TRUE))
{
  if (!is.finite(prevE0)) prevE0 = 100
  resFit <- try(
    nls(
      formula = R_eco ~ fLloydTaylor(RRef, E0, Temp, TRef = TRefFit),
      algorithm = "default", trace = FALSE,
      data = as.data.frame(cbind(R_eco = REco, Temp = temperatureKelvin)),
      start = list(RRef = mean(REco, na.rm = TRUE), E0 = as.vector(prevE0)),
      control = nls.control(maxiter = 20L)
    ),
    silent = TRUE)
  # plot(REco ~ I(temperatureKelvin-273.15) )
  # lines(fLloydTaylor(coef(resFit)["RRef"], coef(resFit)["E0"]
  # , temperatureKelvin, TRef = TRefFit) ~ I(temperatureKelvin-273.15))
  if (inherits(resFit, "try-error")) {
    # stop("debug partGLEstimateTempSensInBounds")
    # plot(REco.V.n	~ temperatureKelvin.V.n)
    E0 <- NA
    sdE0 <- NA
    RRefFit <- NA
  } else {
    E0 <- coef(resFit)["E0"]
    sdE0 <- coef(summary(resFit))["E0", 2]
    RRefFit <- coef(resFit)["RRef"]
  }
  # resFit$convInfo$isConv
  ## details<<
  ## If E0 is out of bounds [50, 400] then report E0 as NA
  if (is.na(E0) || (E0 < 50) || (E0 > 400)) {
    E0 <- NA
    sdE0 <- NA
    RRefFit <- NA
  }
  return(list(E0 = E0, sdE0 = sdE0, TRefFit = TRefFit, RRefFit = RRefFit))
}

#' Estimate Reference temperature from nighttime and given tempsens E0
#' 
#' Estimation of respiration at reference temperature (RRef) and temperature
#' sensitivity (E0) for one window.
#' 
#' @details If there are too few records
#' (n < `controlGLPart$minNRecInDayWindow`) then return NA.
#' 
#' @param dss data.frame with numeric columns NEE, isNight, Temp, Rg
#' @param winInfo one-row data.frame with window information, including iWindow
#' @param prevRes component prevRes from previous result, here not used.
#' @param E0Win data.frame with columns E0 and sdE0, RRefFit, and TRefFit
#' with one row for each window
#' @param controlGLPart list of further default parameters
#' @param TRef numeric scalar of Temperature (degree Celsius) for reference
#' respiration RRef
#' 
#' @seealso [partitionNEEGL()], [partGLEstimateTempSensInBoundsE0Only()]
#' @references 
#' 
#' @return named scalar, RRef
#' @export 
partGLFitNightRespRefOneWindow <- function(
  dss, winInfo, prevRes = list(), E0Win, controlGLPart = partGLControl(), TRef = 15){
  
  isValid <- isValidNightRecord(dss)
  if (sum(isValid) < controlGLPart$minNRecInDayWindow) return(c(RRef = NA_real_))
  
  dssNight <- dss[isValid, ]
  REco <- dssNight$NEE
  E0 <- E0Win$E0[winInfo$iWindow]
  TKref <- 273.15 + TRef # 15degC in Kelvin

  RRef <- if (length(REco) >= 3L) {
    K <- 273.15 + dssNight$Temp
    T_0.n <- 227.13 # -46.02 + 273.15
    TFac <- exp(E0 * (1 / (TKref - T_0.n) - 1 / (K - T_0.n))) # LloydTaylor (1994)
    lm15 <- lm(REco ~ TFac - 1)
    coef(lm15)
  } else {
    fLloydTaylor(E0Win$RRefFit[winInfo$iWindow], E0, TKref,
      TRef = E0Win$TRefFit[winInfo$iWindow])
  }
  c(RRef = max(0, RRef))
}
