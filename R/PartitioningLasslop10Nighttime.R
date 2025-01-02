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
#' @seealso [partGLFitNightTempSensOneWindow()], [partGLSmoothTempSens()], 
#' [partGLFitNightRespRefOneWindow()]
#' 
#' @export 
partGLFitNightTimeTRespSens <- function(
  ds, 
  winSizeRefInDays = 4L, winSizeNight = 12L, 
  winExtendSizes = winSizeNight * c(2L, 4L),
  strideInDays = 2L, 
  isVerbose = TRUE, nRecInDay = 48L, 
  controlGLPart = partGLControl()
) {
  if (isVerbose) 
    message("  Estimating temperature sensitivity from night time NEE ", appendLF = FALSE)
  
  resNight <- simplifyApplyWindows(tmp <- applyWindows(ds,
    partGLFitNightTempSensOneWindow,
    prevRes = data.frame(E0 = NA),
    winSizeRefInDays = winSizeRefInDays,
    winSizeInDays = winSizeNight,
    isVerbose = isVerbose,
    nRecInDay = nRecInDay,
    controlGLPart = controlGLPart
  ))
  iNoSummary <- which(is.na(resNight$E0))
  iExtend <- 1
  
  if (!isTRUE(controlGLPart$isExtendTRefWindow)) winExtendSizes <- numeric(0)
  
  while (length(iNoSummary) && (iExtend <= length(winExtendSizes))) {
    if (isVerbose) 
      message("    increase window size to ", winExtendSizes[iExtend], appendLF = FALSE)
    
    resNightExtend <- simplifyApplyWindows(applyWindows(ds,
      partGLFitNightTempSensOneWindow,
      prevRes = data.frame(E0 = NA),
      winSizeRefInDays = winSizeRefInDays,
      winSizeInDays = winExtendSizes[iExtend],
      isVerbose = isVerbose,
      nRecInDay = nRecInDay,
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
    E0Smooth <- partGLSmoothTempSens(resNight)
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
    partGLFitNightRespRefOneWindow,
    winSizeRefInDays = winSizeRefInDays,
    winSizeInDays = winSizeNight,
    nRecInDay = nRecInDay, 
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

#' getStartRecsOfWindows
#' 
#' compute the starting positions of windows for given size
#' @param nRec numeric scalar: number of records
#' @param winSizeInDays Window size in days
#' @param winSizeRefInDays Window size in days for reference window
#' (e.g. day-Window for night time)
#' @param strideInDays step in days for shifting the window, for alligning
#' usually a factor of winSizeRef
#' @param nRecInDay number of records within one day (for half-hourly data its 48)
#' 
#' @examples getStartRecsOfWindows(1000)
#' @export 
getStartRecsOfWindows <- function(
  nRec, winSizeInDays = winSizeRefInDays, winSizeRefInDays = 4L,
  strideInDays = floor(winSizeRefInDays / 2L), nRecInDay = 48L) {
  
  win_half = winSizeInDays / 2
  nDay <- as.integer(ceiling(nRec / nRecInDay))
  
  # center of the reference window still in records
  nDayLastWindow <- nDay - win_half
  # specifying the day for each record assuming equidistand records
  # iDayOfRec <- ((c(1:nRec)-1L) %/% nRecInDay) + 1L
  # starting days for each reference window
  startDaysRef <- seq(1, nDayLastWindow, strideInDays)
  # assuming equidistant records
  iCentralRec <- 1L + as.integer((startDaysRef - 1L) + win_half) * nRecInDay
  # precomputing the starting and end records for all periods in vectorized way
  # may become negative, negative needed for computation of iRecEnd
  iRecStart0 <- as.integer(iCentralRec - win_half * nRecInDay)
  pmax(1L, iRecStart0) # iRecStart
}

# get_winInfo
#' @export 
get_winInfo <- function(
  nRec, 
  winSizeInDays = winSizeRefInDays, 
  winSizeRefInDays = 4L, strideInDays = floor(winSizeRefInDays / 2L),
  nRecInDay = 48L) 
{
  nDay <- as.integer(ceiling(nRec / nRecInDay))
  nDayLastWindow <- nDay - (winSizeRefInDays / 2)
  
  startDaysRef <- seq(1, nDayLastWindow, strideInDays)
  iCentralRec <- 1L + as.integer((startDaysRef - 1L) + winSizeRefInDays / 2) * nRecInDay

  nWindow <- length(startDaysRef)
  dayStart0 <- as.integer(startDaysRef + winSizeRefInDays / 2 - winSizeInDays / 2)
  
  dayStart <- pmax(1L, dayStart0)
  dayEnd <- pmin(nDay, dayStart0 - 1L + winSizeInDays)

  iRecStart0 <- as.integer(iCentralRec - winSizeInDays / 2 * nRecInDay)
  iRecStart <- pmax(1L, iRecStart0)
  iRecEnd <- pmin(nRec, as.integer(iCentralRec - 1L + winSizeInDays / 2 * nRecInDay))

  data.frame(
    iWindow = 1:nWindow ## << integer: counter of the window
    , dayStart = dayStart ## << integer: starting day of the window
    , dayEnd = dayEnd ## << integer: ending day of the window
    , iRecStart = iRecStart ## << integer: first record number of the window
    , iRecEnd = iRecEnd ## << integer: last record number of the window
    , iCentralRec = iCentralRec ## << integer: central record within the window
    ## assuming equidistant records
  )
}

#' apply a function to several windows of a data.frame
#' 
#' @param ds data.frame to iterate
#' @param FUN function to apply to subsets of the data.frame taking a subset of
#' the data.frame as first argument the second: a one-row data.frame with window
#' information (iWindow, dayStart, dayEnd, iRecStart, iRecEnd, iCentralRec) the
#' third: most recent valid result of calls FUN. Valid is a non-NULL result.
#' @param prevRes initial values for the list that
#' is carried between windows
#' @param isVerbose set to FALSE to suppress messages
#' @param ... further arguments to FUN
#' @inheritParams getStartRecsOfWindows
#' 
#' @export 
applyWindows <- function(
  ds, FUN, prevRes = list(), 
  winSizeInDays = winSizeRefInDays, 
  winSizeRefInDays = 4L, strideInDays = floor(winSizeRefInDays / 2L),
  nRecInDay = 48L, 
  isVerbose = TRUE, ...) 
{
  ## details<<
  # Assumes equidistant rows with nRecInDay records forming one day and
  # reporting full days, i.e. all of the nRecInDay records are in the first day.

  ## details<<
  # In order to have a common reference winSizeRefInDays is given so that
  # results by a different window size correspond to each window of shifting a
  # window of winSizeRefInDays Each window is anchord so that the center equals
  # the center of the reference window.
  # This becomes important when selecting records at the edges of the series.
  info <- get_winInfo(nrow(ds), winSizeInDays, winSizeRefInDays, strideInDays, nRecInDay)  
  nWindow <- nrow(info)
  
  # each will hold a data.frame to be row-bound afterwards (dont know the cols yet)
  res2List <- vector("list", nWindow)
  for (iWindow in 1:nWindow) {
    if (isVerbose) message(", ", info$dayStart[iWindow], appendLF = FALSE)
    startRec <- info$iRecStart[iWindow]
    endRec <- info$iRecEnd[iWindow]
    dsWin <- ds[startRec:endRec, ]
    resFun <- FUN(dsWin, info[iWindow, ], prevRes, ...)
    ## details<< Usually indicate an invalid result by returning NULL.
    ## If one still wants to store results but prevent updating the
    ## \code{prevRes} argument supplied to the next call
    ## then return a list item (or dataframe column) \code{isValid = TRUE}.
    if (length(resFun)) {
      res2List[[iWindow]] <- resFun
      if (!is.list(resFun) || !length(resFun$isValid) || isTRUE(resFun$isValid)) {
        prevRes <- resFun
      }
    }
  }
  if (isVerbose) message("") # LineFeed
  list(winInfo = info, resFUN = res2List)
}

#' simplify the result returned by applyWindows
#' 
#' @details 
#' If FUN returns a named vector or a single-row data.frame,
#' the resFUN result component of applyWindows can be condensed to a data.frame.
#' This result is column-bound to the winInfo result component
#' 
#' @param resApply result of [applyWindows()]
#' @return A single data.frame with columns of winInfo and results of FUN
simplifyApplyWindows <- function(resApply) {
  if (!length(resApply$winInfo)) {
    return(resApply$winInfo)
  }
  ansFUN <- if (is.data.frame(resApply$resFUN[[1]])) {
    bind_rows(resApply$resFUN)
  } else {
    do.call(rbind, resApply$resFUN)
  }
  cbind(resApply$winInfo, ansFUN)
}

#' Smoothes time development of E0
#' @param E0Win data.frame with columns E0 and sdE0, RRefFit, and TRefFit
#' with one row for each window
partGLSmoothTempSens <- function(E0Win) {
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
