# ' @param LRCFitConvergenceTolerance convergence criterion for rectangular light
# ' response curve fit. If relative improvement of reducing residual sum of
# ' squares between predictions and observations is less than this criterion,
# ' assume convergence. Decrease to get more precise parameter estimates,
# ' Increase for speedup.

#' @export
partGLControl <- function(
    ### Default list of parameters for Lasslop 2010 daytime flux partitioning
    LRCFitConvergenceTolerance = 1e-3 ## << convergence criterion for rectangular
    ## light response curve fit.
    ## If relative improvement of reducing residual sum of squares between
    ## predictions and
    ## observations is less than this criterion, assume convergence.
    ## Decrease to get more precise parameter estimates, Increase for speedup.
    , nLRCFitConvergenceTolerance = 1e-3 ## << convergence criterion for
    ## nonrectangular light response curve fit.
    ## Here its a factor of machine tolerance.
    , nBootUncertainty = 30L ## << number of bootstrap samples for
    ## estimating uncertainty.
    ## Set to zero to derive uncertainty from curvature of a single fit
    , minNRecInDayWindow = 10L ## << Minimum number of data points
    ## for regression
    , isAssociateParmsToMeanOfValids = TRUE ## << set to FALSE to
    ## associate parameters to
    ## the first record of the window for interpolation
    ## instead of mean across valid records inside a window
    , isLasslopPriorsApplied = TRUE ## << set to TRUE to apply strong fixed
    ## priors on LRC fitting.
    ## Returned parameter estimates claimed valid for some case where not
    ## enough data was available
    , isUsingLasslopQualityConstraints = FALSE ## << set to TRUE to avoid
    ## quality constraints additional to Lasslop 2010
    , isSdPredComputed = TRUE ## << set to FALSE to avoid computing
    ## standard errors
    ## of Reco and GPP for small performance increase
    , isFilterMeteoQualityFlag = FALSE ## << set to TRUE to use only records
    ## where quality flag
    ## of meteo drivers (radiation, temperature, VPD) is zero, i.e.
    ## non-gapfilled for parameter estimation.
    ## For prediction, the gap-filled value is used always, to produce
    ## predictions also for gaps.
    , isBoundLowerNEEUncertainty = TRUE ## << set to FALSE to avoid adjustment
    ## of very low uncertainties before
    ## day-Time fitting that avoids the high leverage those records with
    ## unreasonable low uncertainty.
    , fixedTRefAtNightTime = NA ## << if a finite value (degree Centigrade)
    ## is given, it is used instead of median data temperature as reference
    ## temperature in estimation of temperature sensitivity from night data
    , isExtendTRefWindow = TRUE ## << set to FALSE to avoid successively
    ## extending the night-time window
    ## in order to estimate a temperature sensitivity where previous estimates
    ## failed
    , smoothTempSensEstimateAcrossTime = TRUE ## << set to FALSE to use
    ## independent estimates of temperature
    ## sensitivity on each windows instead of a vector of E0 that is
    ## smoothed over time
    , isNeglectPotRadForNight = FALSE ## << set to TRUE to not use potential
    ## radiation in determining night-time data.
    , NRHRfunction = FALSE ## << deprecated: Flag if TRUE use the NRHRF
    ## for partitioning; Now use \code{lrcFitter = NonrectangularLRCFitter()}
    , isNeglectVPDEffect = FALSE ## << set to TRUE to avoid using VPD in the
    ## computations. This may help when VPD is rarely measured.
    , isRefitMissingVPDWithNeglectVPDEffect = TRUE ## << set to FALSE to avoid
    ## repeating estimation
    ## with \code{isNeglectVPDEffect = TRUE} trying to predict when VPD
    ## is missing
    , fixedTempSens = data.frame( ## << data.frame
      ## of one row or nRow = nWindow
      ## corresponding to return value of \code{partGL_FitNight_E0_RRef}
      ## While column \code{RRef} is used only as a  prior and initial value for
      ## the daytime-fitting and can be NA,
      ## \code{E0} is used as given temperature sensitivity and varied according
      ## to \code{sdE0} in the bootstrap.
      E0 = NA_real_, sdE0 = NA_real_, RRef = NA_real_
    ),
    replaceMissingSdNEEParms = c(perc = 0.2, minSd = 0.7) ## << parameters for
    ## replacing missing standard deviation of NEE.
    ## see \code{replaceMissingSdByPercentage}.
    ## Default sets missing uncertainty to 20% of NEE but at least 0.7
    ## flux-units (usually mumol CO2 / m2 / s).
    ## Specify c(NA, NA) to avoid replacing missings in standard deviation of
    ## NEE and to omit those records from LRC fit.
    , neglectNEEUncertaintyOnMissing = FALSE ## << If set to TRUE: if there are
    ## records with missing uncertainty of NEE inside one window,
    ## set all uncertainties to 1.
    ## This overrules option replaceMissingSdNEEParms.
    , minPropSaturation = NA ## << quality criterion for sufficient data
    ## in window. If GPP prediction of highest PAR of window is less than
    ## minPropSaturation * (GPP at light-saturation, i.e. beta)
    ## this indicates that PAR is not sufficiently high to constrain the
    ## shape of the LRC
    , useNightimeBasalRespiration = FALSE ## << set to TRUE to estimate
    ## nighttime respiration based on basal respiration estimated on
    ## nighttime data instead of basal respiration estimated from daytime
    ## data. This implements the modified daytime method from
    ## Keenan 2019 (doi:10.1038/s41559-019-0809-2)
    ) {
  ## author<< TW
  ## seealso<< \code{\link{partitionNEEGL}}
  ## description<<
  ## For highest compatibility to the pvWave code of G.Lasslop
  ## (used by first BGC-online tool)
  ## see function \code{\link{partGLControlLasslopCompatible}}.
  if (NRHRfunction) {
    stop(
      "option 'NRHRfunction' is deprecated.", " Use instead in partitionNEEGL argument:",
      "lrcFitter = NonrectangularLRCFitter()"
    )
  }
  if (isTRUE(neglectNEEUncertaintyOnMissing)) {
    replaceMissingSdNEEParms <- c(NA, NA)
  }

  listk(
    LRCFitConvergenceTolerance,
    nLRCFitConvergenceTolerance,
    nBootUncertainty,
    minNRecInDayWindow,
    isAssociateParmsToMeanOfValids,
    isLasslopPriorsApplied,
    isUsingLasslopQualityConstraints,
    isSdPredComputed,
    isFilterMeteoQualityFlag,
    isBoundLowerNEEUncertainty,
    fixedTRefAtNightTime,
    isExtendTRefWindow,
    smoothTempSensEstimateAcrossTime,
    isNeglectPotRadForNight,
    isNeglectVPDEffect,
    isRefitMissingVPDWithNeglectVPDEffect,
    fixedTempSens,
    replaceMissingSdNEEParms,
    neglectNEEUncertaintyOnMissing,
    minPropSaturation,
    useNightimeBasalRespiration
  )
}
attr(partGLControl, "ex") <- function() {
  partGLControl(nBootUncertainty = 40L)
}

#' @export
partGLControlLasslopCompatible <- function(
    ### Daytime flux partitioning parms compatible with with the pvWave
    nBootUncertainty = 0L ## << 0: Derive uncertainty from
    ## curvature of a single fit, neglecting the uncertainty of previously
    ## estimated temperature sensitivity, E0
    , minNRecInDayWindow = 10L ## << Minimum number of 10 valid records
    ## for regression in a single window
    , isAssociateParmsToMeanOfValids = FALSE ## << associate parameters to
    ## the first record of the window for interpolation instead of mean across
    ## valid records inside a window
    , isLasslopPriorsApplied = TRUE ## << Apply fixed Lasslop priors
    ## in LRC fitting.
    , isUsingLasslopQualityConstraints = TRUE ## << avoid quality constraints
    ## additional to the ones in Lasslop 2010
    , isBoundLowerNEEUncertainty = FALSE ## << FALSE: avoid adjustment of very
    ## low uncertainties before
    ## day-Time fitting that avoids the high leverage those records with
    ## unreasonable low uncertainty.
    , fixedTRefAtNightTime = 15 ## << use fixed (degree Centigrade)
    ## temperature sensitivity
    ## instead of median data temperature as reference temperature in
    ## estimation of temperature sensitivity from night data
    , isExtendTRefWindow = FALSE ## << avoid successively extending the
    ## night-time window
    ## in order to estimate a temperature sensitivity where previous
    ## estimates failed
    , smoothTempSensEstimateAcrossTime = FALSE ## << FALSE: use independent
    ## estimates of temperature
    ## sensitivity on each windows instead of a vector of E0 that is
    ## smoothed over time
    , isRefitMissingVPDWithNeglectVPDEffect = FALSE ## << FALSE: avoid
    ## repeating estimation with \code{isNeglectVPDEffect = TRUE}
    , minPropSaturation = NA ## << NA: avoid quality constraint of
    ## sufficient saturation in data
    ## This option is overruled, i.e. not considered, if option
    ## isUsingLasslopQualityConstraints = TRUE.
    , isNeglectVPDEffect = FALSE ## << FALSE: do not neglect VPD effect
    , replaceMissingSdNEEParms = c(NA, NA) ## << do not replace missing NEE,
    ## but see option
    , neglectNEEUncertaintyOnMissing = TRUE ## << if there are records with
    ## missing uncertainty of NEE inside one window, set all sdNEE to 1.
    ## This overrules option replaceMissingSdNEEParms.
    , ... ## << further arguments to \code{\link{partGLControl}}
    ) {
  ## seealso<< \code{\link{partGLControl}}
  partGLControl(
    nBootUncertainty = nBootUncertainty,
    minNRecInDayWindow = minNRecInDayWindow,
    isAssociateParmsToMeanOfValids = isAssociateParmsToMeanOfValids,
    isLasslopPriorsApplied = isLasslopPriorsApplied,
    isUsingLasslopQualityConstraints = isUsingLasslopQualityConstraints,
    isBoundLowerNEEUncertainty = isBoundLowerNEEUncertainty,
    fixedTRefAtNightTime = fixedTRefAtNightTime,
    isExtendTRefWindow = isExtendTRefWindow,
    smoothTempSensEstimateAcrossTime = smoothTempSensEstimateAcrossTime,
    isNeglectVPDEffect = isNeglectVPDEffect,
    isRefitMissingVPDWithNeglectVPDEffect = isRefitMissingVPDWithNeglectVPDEffect,
    replaceMissingSdNEEParms = replaceMissingSdNEEParms,
    neglectNEEUncertaintyOnMissing = neglectNEEUncertaintyOnMissing,
    minPropSaturation = minPropSaturation,
    ...
  )
}
attr(partGLControlLasslopCompatible, "ex") <- function() {
  partGLControlLasslopCompatible()
}
