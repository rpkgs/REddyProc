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
  win_half <- winSizeInDays / 2
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
    nRecInDay = 48L) {
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
    ds, FUN, prevRes = list(), winInfo = get_winInfo(nrow(ds)),
    isVerbose = TRUE, ...) {
  ## details<<
  # Assumes equidistant rows with nRecInDay records forming one day and
  # reporting full days, i.e. all of the nRecInDay records are in the first day.

  ## details<<
  # In order to have a common reference winSizeRefInDays is given so that
  # results by a different window size correspond to each window of shifting a
  # window of winSizeRefInDays Each window is anchord so that the center equals
  # the center of the reference window.
  # This becomes important when selecting records at the edges of the series.

  # winInfo <- get_winInfo(nrow(ds), winSizeInDays, winSizeRefInDays, strideInDays, nRecInDay)
  nWindow <- nrow(winInfo)

  # each will hold a data.frame to be row-bound afterwards (dont know the cols yet)
  res2List <- vector("list", nWindow)
  for (iWindow in 1:nWindow) {
    if (isVerbose) message(", ", winInfo$dayStart[iWindow], appendLF = FALSE)
    startRec <- winInfo$iRecStart[iWindow]
    endRec <- winInfo$iRecEnd[iWindow]
    dsWin <- ds[startRec:endRec, ]
    resFun <- FUN(dsWin, winInfo[iWindow, ], prevRes, ...)
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
  list(winInfo = winInfo, resFUN = res2List)
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
