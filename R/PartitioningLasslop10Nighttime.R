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


#' Partitioning of NEE into GPP and RE using the Lasslop et al. (2010) method
#' 
#' Estimate temperature sensitivity parameters for successive periods
#' 
#' @details 
#' Window sizes are successively increased to `winExtendSizes` in order to obtain
#' parameter estimates where fits with smaller window size failed. This behaviour
#' of repeated attempts can be avoid by setting `controlGLPart$isExtendTRefWindow = FALSE`.
#' 
#' @param ds data.frame with numeric columns NEE, sdNEE, Temp (degC),
#' VPD, Rg, and logical columns isNight and isDay
#' @param winSizeRefInDays Window size in days for daytime for referencing
#' windows to the same base
#' @param winSizeNight Window size in days for nighttime fits
#' @param winExtendSizes successively increased nighttime windows, to obtain a
#' night-time fit
#' @param strideInDays step in days for shifting the windows
#' @param nRecInDay number of records within one day (for half-hourly data its 48)
#' @param isVerbose set to FALSE to suppress messages
#' @param controlGLPart list of further default parameters
#' 
#' @seealso [partGL_FitNight_1win_E0()], [partGL_smoothE0()], 
#' [partGL_FitNight_1win_RRef()]
#' 
#' @export 
partGL_FitNight_E0_RRef <- function(
  ds, 
  winSizeRefInDays = 4L, winSizeNight = 12L, 
  winExtendSizes = winSizeNight * c(2L, 4L),
  strideInDays = 2L, nRecInDay = 48L, 
  isVerbose = TRUE, 
  controlGLPart = partGLControl()
) {
  if (isVerbose) 
    message("  Estimating temperature sensitivity from night time NEE ", appendLF = FALSE)
  
  winInfo <- get_winInfo(nrow(ds),
    winSizeRefInDays = winSizeRefInDays,
    winSizeInDays = winSizeNight,
    nRecInDay = nRecInDay
  )
  
  resNight <- simplifyApplyWindows(tmp <- applyWindows(ds,
    partGL_FitNight_1win_E0,
    prevRes = data.frame(E0 = NA), winInfo = winInfo,
    isVerbose = isVerbose, controlGLPart = controlGLPart))
  iNoSummary <- which(is.na(resNight$E0))
  iExtend <- 1
  
  if (!isTRUE(controlGLPart$isExtendTRefWindow)) winExtendSizes <- numeric(0)
  
  while (length(iNoSummary) && (iExtend <= length(winExtendSizes))) {
    if (isVerbose) 
      message("    increase window size to ", winExtendSizes[iExtend], appendLF = FALSE)
    
    .winInfo <- get_winInfo(nrow(ds),
      winSizeRefInDays = winSizeRefInDays,
      winSizeInDays = winExtendSizes[iExtend],
      nRecInDay = nRecInDay
    )
    resNightExtend <- simplifyApplyWindows(applyWindows(ds,
      partGL_FitNight_1win_E0,
      prevRes = data.frame(E0 = NA), winInfo = .winInfo,
      isVerbose = isVerbose,
      controlGLPart = controlGLPart
    )) # T_ref固定(median(temp)), 优化得到`E0`, `RRef`
    resNight[iNoSummary, ] <- resNightExtend[iNoSummary, ]
    # all(resNight$iCentralRec == resNightExtend$iCentralRec)
    iNoSummary <- which(is.na(resNight$E0))
    iExtend <- iExtend + 1L
  }
  
  # remember E0 and sdE0 before overidden by smoothing
  nFiniteE0 <- sum(is.finite(resNight$E0))
  if ((nFiniteE0 < 5) && (nFiniteE0 < 0.1 * nrow(resNight))) {
    stop("Estimated valid temperature sensitivity for only ", nFiniteE0,
      " windows. Stopping.")
  }
  resNight$E0Fit <- resNight$E0
  resNight$sdE0Fit <- resNight$sdE0
  E0Smooth <- if (isTRUE(controlGLPart$smoothTempSensEstimateAcrossTime)) {
    if (isVerbose) message("  Smoothing temperature sensitivity estimates")
    E0Smooth <- partGL_smoothE0(resNight)
  } else {
    E0Smooth <- resNight
    iNonFiniteE0 <- which(!is.finite(E0Smooth$E0))
    # set uncertainty to the 90% quantile of the distribution of uncertainties
    E0Smooth$sdE0[iNonFiniteE0] <- quantile(E0Smooth$sdE0, 0.9, na.rm = TRUE)
    # fill NA with value from previous window
    E0Smooth$E0 <- fillNAForward(E0Smooth$E0)
    E0Smooth
  }
  # now all E0 and sdE0 are defined

  if (isVerbose) message(
      "  Estimating respiration at reference temperature for smoothed temperature",
      " sensitivity from night time NEE ", appendLF = FALSE)
  
  resRef15 <- simplifyApplyWindows(tmp <- applyWindows(ds,
    partGL_FitNight_1win_RRef, winInfo = winInfo, 
    E0Win = E0Smooth,
    controlGLPart = controlGLPart,
    isVerbose = isVerbose
  )) # E0固定，T_ref固定(15℃)，得到平滑的`RRef`
  # fill those NA caused by not enough night-time records
  E0Smooth$RRef <- fillNAForward(resRef15$RRef,
    firstValue = E0Smooth$RRef[which(is.finite(E0Smooth$RRef))[1]])
  # on NA at the start of the series, take the first  finite value
  
  ## value<< data.frame with columns of \code{winInfo}
  ## from applyWindows, E0, sdE0, RRef
  E0Smooth
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
#' @seealso [partitionNEEGL()], [partGL_FitNight_1win_E0_kernel()]
partGL_FitNight_1win_E0 <- function(
    dss, winInfo, prevRes, isVerbose = TRUE, nRecInDay = 48L,
    controlGLPart = partGLControl()) {
  isValid <- isValidNightRecord(dss)
  # check that there are enough night and enough day-values for fitting,
  # else continue with next window
  if (sum(isValid) < controlGLPart$minNRecInDayWindow) {
    return(data.frame(E0 = NA_real_, sdE0 = NA_real_, TRefFit = NA_real_, RRefFit = NA_real_))
  }

  dssNight <- dss[isValid, ]

  TRefFitCelsius <- if (length(controlGLPart$fixedTRefAtNightTime) &&
    is.finite(controlGLPart$fixedTRefAtNightTime)
  ) {
    controlGLPart$fixedTRefAtNightTime
  } else {
    median(dssNight$Temp, na.rm = TRUE)
  }
  resNightFit <- partGL_FitNight_1win_E0_kernel(
    dssNight$NEE,
    fConvertCtoK(dssNight$Temp),
    prevE0 = prevRes$E0, TRefFit = fConvertCtoK(TRefFitCelsius)
  )
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
partGL_FitNight_1win_E0_kernel <- function(
    REco, temperatureKelvin, prevE0 = NA,
    TRefFit = median(temperatureKelvin, na.rm = TRUE)) {
  #
  if (!is.finite(prevE0)) prevE0 <- 100
  resFit <- try(
    nls(
      formula = R_eco ~ fLloydTaylor(RRef, E0, Temp, TRef = TRefFit),
      algorithm = "default", trace = FALSE,
      data = as.data.frame(cbind(R_eco = REco, Temp = temperatureKelvin)),
      start = list(RRef = mean(REco, na.rm = TRUE), E0 = as.vector(prevE0)),
      control = nls.control(maxiter = 20L)
    ),
    silent = TRUE
  )
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
#' @seealso [partitionNEEGL()], [partGL_FitNight_1win_E0_kernel()]
#'
#' @return named scalar, RRef
#' @export
partGL_FitNight_1win_RRef <- function(
    dss, winInfo, prevRes = list(), E0Win, controlGLPart = partGLControl(), TRef = 15) {
  isValid <- isValidNightRecord(dss)
  if (sum(isValid) < controlGLPart$minNRecInDayWindow) {
    return(c(RRef = NA_real_))
  }

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
      TRef = E0Win$TRefFit[winInfo$iWindow]
    )
  }
  c(RRef = max(0, RRef))
}

#' Smoothes time development of E0
#' @param E0Win data.frame with columns E0 and sdE0, RRefFit, and TRefFit
#' with one row for each window
partGL_smoothE0 <- function(E0Win) {
  # E0Win$E0[1] <- NA
  # TODO return NA in the first place where the previous window was used
  E0Win$E0[c(FALSE, (diff(E0Win$E0) == 0))] <- NA
  # mlegp does not work for long time series (Trevors report with Havard data)
  # hence, smooth one year at a time
  E0Win$year <- ceiling(E0Win$dayStart / 365)
  yr <- E0Win$year[1]

  E0_smooth <- function(yr) {
    E0WinYr <- subset(E0Win, year == yr)
    isFiniteE0 <- is.finite(E0WinYr$E0)
    E0WinFinite <- E0WinYr[isFiniteE0, ]
    if (!nrow(E0WinFinite)) {
      warning(
        "No respiration-temperature relationship for any period of the ", yr, "th year. ",
        "Using mean temperature sensitivity of other years for this year.")
    } else {
      if (!requireNamespace("mlegp")) {
        stop(
          "package mlegp is required for Lasslop daytime partitioning.",
          "Unfortunately, it is removed from CRAN.",
          "In order to still use daytime partitioning, please, install mlegp ",
          "(maybe an archived version) from ",
          "https://cran.r-project.org/package=mlegp")
      }
      output <- capture.output(
        gpFit <- mlegp::mlegp(
          X = E0WinFinite$iCentralRec, Z = E0WinFinite$E0, 
          nugget = E0WinFinite$sdE0^2)
        # gpFit <- mlegp(X = E0WinFinite$iCentralRec, Z = E0WinFinite$E0
        #              , nugget = (E0WinFinite$sdE0 * 2)^2, nugget.known = 1L)
      )
      pred1 <- predict(gpFit, matrix(E0WinYr$iCentralRec, ncol = 1), se.fit = TRUE)
      nuggetNewObs <- quantile(gpFit$nugget, 0.9)
      nugget <- rep(nuggetNewObs, nrow(E0WinYr))
      nugget[isFiniteE0] <- gpFit$nugget
      E0WinYr$E0 <- as.vector(pred1$fit)
      E0WinYr$sdE0 <- as.vector(pred1$se.fit) + unname(sqrt(nugget))
    }
    E0WinYr
  }

  E0WinYrs <- lapply(unique(E0Win$year), E0_smooth) %>% do.call(rbind, .)
  ## value<< dataframe E0Win with updated columns E0 and sdE0
  
  iNoFiniteE0 <- which(!is.finite(E0WinYrs$E0))
  if (length(iNoFiniteE0)) {
    if (length(iNoFiniteE0) == nrow(E0WinYrs)) {
      stop("Could not estimate respiration~temperature relationship.")
    }
    E0WinYrs$E0[iNoFiniteE0] <- mean(E0WinYrs$E0[-iNoFiniteE0])
    E0WinYrs$sdE0[iNoFiniteE0] <- quantile(E0WinYrs$sdE0[-iNoFiniteE0], 0.9) * 1.5
  }
  E0WinYrs
}

# .tmp.f <- function() {
#   # from recover at fitting E0
#   temp <- temperatureKelvin - 273.15
#   plot(REco ~ temp)
#   lines(fLloydTaylor(RRefFit, E0, temperatureKelvin, TRef = TRefFit) ~ temp)
# }

# .tmp.f <- function() {
#   plot(E0 ~ iCentralRec, E0Win, type = "l")
#   arrows(E0WinFinite$iCentralRec, E0WinFinite$E0 - E0WinFinite$sdE0,
#     y1 = E0WinFinite$E0 + E0WinFinite$sdE0, length = 0, col = "grey"
#   )
#   plot(E0 ~ iCentralRec, E0WinYrs, type = "l")
#   plot(sdE0 ~ iCentralRec, E0WinYrs, type = "l")
#   #
#   E0Win$day <- (E0Win$iCentralRec - 1) / 48 + 1
#   E0WinFinite$day <- (E0WinFinite$iCentralRec - 1) / 48 + 1
#   plot(E0WinFinite$E0 ~ E0WinFinite$day)
#   plot(E0Win$E0Fit ~ E0Win$day)
#   points(E0Win$E0 ~ E0Win$day, col = "blue", type = "b", lty = "dotted")
#   # points(E0Win$E0 ~ E0Win$day, col = "blue")
#   # arrows(E0Win$day, E0Win$E0Fit, y1 = E0Win$E0, col = "grey", length = 0.1)
#   lines(I(E0Win$E0 + 1.06 * E0Win$sdE0) ~ E0Win$day, col = "lightblue")
#   lines(I(E0Win$E0 - 1.06 * E0Win$sdE0) ~ E0Win$day, col = "lightblue")
#   #
#   E0Win$E0 <- E0Win$E0Fit
#   E0Win$sdE0 <- E0Win$sdE0Fit
#   E0Win$E0[c(FALSE, (diff(E0Win$E0) == 0))] <- NA # TODO ret NA in the first place
#   #
#   E0Win$isFiniteE0 <- isFiniteE0
#   E0Win[E0Win$E0 < 80, ]
#   E0Win[65:75, ]
#   subset(E0WinFinite, iWindow %in% 65:75)

#   E0Win$sdE0 / E0Win$E0
# }
