#' Partition NEE fluxes into GPP and Reco using the daytime method.
#'
#' @details The daytime-based partitioning of measured net ecosystem fluxes into
#' gross primary production (GPP) and ecosystem respiration (Reco).
#'
#' The fit to the light-response-curve is done by default using the Rectangular
#' hyperbolic function, as in Lasslop et al. (2010).
#'
#' Alternative fittings can be used by providing the corresponding subclass of
#' `LightResponseCurveFitter-class` to `lrcFitter` argument.
#'
#' see `LRC_predictGPP()`
#'
#' @param ds dataset with all the specified input columns and full days in equidistant times
#' @param NEEVar Variable of NEE
#' @param TempVar Filled air or soil temperature variable (degC)
#' @param VPDVar Filled Vapor Pressure Deficit - VPD - (hPa)
#' @param RadVar Filled radiation variable
#' @param suffix string inserted into column names before identifier for NEE
#' column defaults
#' @param ... further arguments to [partGLExtractStandardData()], such as `PotRadVar`
#' @param controlGLPart further default parameters, see [partGLControl()]
#' @param isVerbose set to FALSE to suppress output messages
#' @param nRecInDay number of records within one day (for half-hourly data its 48)
#' @param lrcFitter R5 class instance responsible for fitting the light response curve.
#' Current possibilities are `RectangularLRCFitter()`, `NonrectangularLRCFitter()`,
#' and `LogisticSigmoidLRCFitter()`.
#'
#' @references
#' 1. Lasslop G, Reichstein M, Papale D, et al. (2010) Separation of net
#'    ecosystem exchange into assimilation and respiration using a light
#'    response curve approach: critical issues and global evaluation. Global
#'    Change Biology, Volume 16, Issue 1, Pages 187-208
#'
#' @return A data.frame with columns:
#' - `Reco_DT_<suffix>`: predicted ecosystem respiration: mumol CO2/m2/s
#' - `GPP_DT_<suffix>`: predicted gross primary production mumol CO2/m2/s
#'
#' - `<LRC>`: Further light response curve (LRC) parameters and their standard
#' deviation depend on the used LRC (e.g. for the non-rectangular LRC see
#' `NonrectangularLRCFitter_getParamNames`). They are estimated for windows and
#' are reported with the first record of the window.
#' 
#' - `FP_VARnight`   : NEE filtered for nighttime records (others NA)
#' - `FP_VARday`     : NEE filtered for daytime records (others NA)
#' - `NEW_FP_Temp`   : temperature after filtering for quality flag, degree Celsius
#' - `NEW_FP_VPD`    : vapour pressure deficit after filtering for quality flag, hPa
#' - `FP_RRef_Night` : basal respiration estimated from nighttime (W / m2)
#' - `FP_qc`         : quality flag:
#'     + 0 : good parameter fit,
#'     + 1 : some parameters out of range, required refit,
#'     + 2 : next parameter estimate is more than two weeks away
#' - `FP_dRecPar`    : records until or after closest record that has a parameter estimate associated
#' - `FP_errorcode`  : information why LRC-fit was not successful or was rejected, see result of [LRC_fitLRC()]
#' - `FP_GPP2000`    : predicted GPP at VPD = 0 and PAR = 2000 : a surrogate for maximum photosynthetic capacity
#' - `FP_OPT_VPD`    : list object of fitting results including `iOpt` and `covParms`
#' - `FP_OPT_NoVPD`  : list object of fitting results with fit neglecting the VPD effect
#' 
#' @seealso [partGL_FitNight_E0_RRef()]
#' @export
partitionNEEGL <- function(
    ds,
    NEEVar = if (!missing(NEEVar.s)) NEEVar.s else paste0("NEE", suffixDash, "_f"),
    TempVar = if (!missing(TempVar.s)) TempVar.s else "Tair_f",
    VPDVar = if (!missing(VPDVar.s)) VPDVar.s else "VPD_f",
    RadVar = if (!missing(RadVar.s)) RadVar.s else "Rg_f",
    suffix = if (!missing(Suffix.s)) Suffix.s else "",
    NEEVar.s, TempVar.s, VPDVar.s, RadVar.s, Suffix.s # deprecated
    , ..., 
    controlGLPart = partGLControl(), isVerbose = TRUE, nRecInDay = 48L,
    lrcFitter = RectangularLRCFitter()) 
{
  varNamesDepr <- c("NEEVar.s", "TempVar.s", "VPDVar.s", "RadVar.s", "Suffix.s")
  varNamesNew <- c("NEEVar", "TempVar", "VPDVar", "RadVar", "suffix")
  iDepr <- which(!c(
    missing(NEEVar.s), missing(TempVar.s), missing(VPDVar.s),
    missing(RadVar.s), missing(Suffix.s)
  ))
  if (length(iDepr)) 
    warning("Argument names ", varNamesDepr[iDepr], " have been deprecated.",
      " Please, use instead ", varNamesNew[iDepr])
  
  suffixDash <- paste0((if (fCheckValString(suffix)) "_" else ""), suffix) # used to compute default NEEVar
  if (isVerbose) {
    message("Start daytime flux partitioning for variable ",
      NEEVar, " with temperature ", TempVar, ".")
  }
  dsR <- partGLExtractStandardData(ds,
    NEEVar = NEEVar, TempVar = TempVar,
    VPDVar = VPDVar, RadVar = RadVar, suffix = suffix,
    ..., controlGLPart = controlGLPart)

  dsAns0 <- tibble(
    FP_VARnight     = rep(NA_real_, nrow(ds)) # NEE filtered for nighttime records (others NA)
    , FP_VARday     = NA_real_                # NEE filtered for daytime records (others NA)
    , NEW_FP_Temp   = dsR$Temp                # temperature after filtering for quality flag, degree Celsius
    , NEW_FP_VPD    = dsR$VPD                 # vapour pressure deficit after filtering for quality flag, hPa
    , FP_RRef_Night = NA_real_                # basal respiration estimated from nighttime (W / m2)
    , FP_qc         = NA_integer_             
    , FP_dRecPar    = NA_integer_             
    , FP_errorcode  = NA_integer_             # information why LRC-fit was not successful or was rejected
    , FP_GPP2000    = NA_real_                # predicted GPP at VPD = 0 and PAR = 2000
    , FP_OPT_VPD    = vector(mode = "list", length = nrow(ds)) 
    , FP_OPT_NoVPD  = vector(mode = "list", length = nrow(ds))
  )

  # append LRC parameter result columns
  lrcParNames <- lrcFitter$getParamNames()
  lrcParNames <- c(lrcParNames, paste0(lrcParNames, "_sd"))
  FP_lrcParNames <- paste0("FP_", lrcParNames)
  tmp <- matrix(NA_real_, nrow = nrow(ds), ncol = length(FP_lrcParNames),
    dimnames = list(NULL, FP_lrcParNames))

  dsAns <- cbind(dsAns0, tmp)

  .units = attr(dsR$NEE, "units")
  .varnames = attr(dsR$NEE, "varnames")
  dsAns$FP_VARnight <- ifelse(dsR$isNight, dsR$NEE, NA) %>% 
    set_attributes(list(varnames = paste0(.varnames, "_night"), units = .units))
  dsAns$FP_VARday <- ifelse(dsR$isDay, dsR$NEE, NA) %>% 
    set_attributes(list(varnames = paste0(.varnames, "_day"), units = .units))
  # Estimate Parameters of light response curve: R_ref, alpha, beta and k
  # according to Table A1 (Lasslop et al., 2010)

  isUsingFixedTempSens <- length(controlGLPart$fixedTempSens) &&
    all(is.finite(controlGLPart$fixedTempSens$E0)) &&
    all(is.finite(controlGLPart$fixedTempSens$sdE0))
  
  dsTempSens <- if (!isUsingFixedTempSens) {
    partGL_FitNight_E0_RRef(dsR,
      nRecInDay = nRecInDay,
      controlGLPart = controlGLPart,
      isVerbose = isVerbose)
  } else {
    dsTempSens <- controlGLPart$fixedTempSens
    nWindow <- length(getStartRecsOfWindows(nrow(dsR), nRecInDay = nRecInDay))
    # if one row is given, use it for all windows
    if (nrow(dsTempSens) == 1L) {
      dsTempSens <- do.call(rbind,
        lapply(1:nWindow, function(i) dsTempSens))
    }
    if (nrow(dsTempSens) != nWindow) {
      stop(
        "when providing controlGLPart$fixedTempSens then it must be 1 row",
        "or rows must macht the number of windows (", nWindow, "), but was ",
        nrow(dsTempSens))
    }
    # if no RRef is given, estimate from nighttime data
    if (!length(dsTempSens$RRef) || all(is.na(dsTempSens$RRef))) {
      resRef15 <- simplifyApplyWindows(applyWindows(dsR,
        partGL_FitNight_1win_RRef,
        nRecInDay = nRecInDay,
        E0Win = dsTempSens,
        controlGLPart = controlGLPart,
        isVerbose = isVerbose
      ))
      dsTempSens$RRef <- resRef15$RRef
    }
    
    dsTempSens$RRef <- fillNAForward(dsTempSens$RRef,
      firstValue = dsTempSens$RRef[which(is.finite(dsTempSens$RRef))[1]])
    dsTempSens
  }
  ## seealso<< \code{partGLFitLRCWindows}
  resParms <- resParmsWithVPD <- partGLFitLRCWindows(dsR,
    nRecInDay = nRecInDay,
    dsTempSens = dsTempSens,
    controlGLPart = controlGLPart,
    lrcFitter = lrcFitter,
    isVerbose = isVerbose
  )
  # if no windows was fitted, return error. Else error on missing columns
  # (no parameter columns returned)
  if (sum(resParms$convergence == 0) == 0) {
    stop("could not fit a single window in dayTime partitioning. Check the data.")
  }
  # append parameter fits to the central record of day window
  # iGood <- which(resLRC$summary$parms_out_range == 0L)
  # resLRC provides parameters only for a subset of rows. For each row in the
  # original data.frame
  # The parameter estimates are associated for interpolation of predicted fluxes
  # either to the central record of the window,
  # or the record with the time corresponding to the mean of all valid records
  # in the window
  # default is isAssociateParmsToMeanOfValids = TRUE
  # (double check partGLControl argument)
  colNameAssoc <- ifelse(isTRUE(controlGLPart$isAssociateParmsToMeanOfValids), 
    "iMeanRec", "iCentralRec")
  # for the output, always report at central record
  dsAns[resParms$iCentralRec, c(
    "FP_RRef_Night", "FP_qc", "FP_errorcode",
    FP_lrcParNames, "FP_GPP2000", "FP_OPT_VPD"
  )] <-
    resParms[, c("RRef_night", "parms_out_range", "convergence", lrcParNames, "GPP2000", "resOpt")]
  ## seealso<< \code{partGLInterpolateFluxes}
  dsAnsFluxes <- partGLInterpolateFluxes(
    ds[[RadVar]], ds[[VPDVar]],
    ds[[TempVar]] # prediction using non-filtered, i.e. gap-filled
    , resParms,
    controlGLPart = controlGLPart,
    lrcFitter = lrcFitter,
    isVerbose = isVerbose
  )
  if (!controlGLPart$isNeglectVPDEffect &&
    controlGLPart$isRefitMissingVPDWithNeglectVPDEffect) {
    ## details<<
    ## While the extrapolation uses filled data, the parameter optimization
    ## may use only measured data, i.e. with specified quality flag.
    ## Even with using filled VPD, there may be large gaps that have not been
    ## filled.
    ## With the common case where VPD is missing for fitting the LRC, by default
    ## (with \code{controlGLPart$isRefitMissingVPDWithNeglectVPDEffect = TRUE})
    ## is to redo the estimation of LRC parameters with neglecting the VPD-effect.
    ## Next, in the predictions (rows) with missing VPD are then replaced
    ## with predictions
    ## based on LRC-fits that neglected the VPD effect.
    iNAVPD <- which(is.na(ds[[VPDVar]] & is.na(dsAnsFluxes$GPP)))
    if (length(iNAVPD)) {
      message("  could not predict GPP in ", length(iNAVPD),
        " cases due to missing VPD.")
      message("    Therefore refitting LightResponseCurve with option ",
        "isNeglectVPDEffect = TRUE")
      ctrlNeglectVPD <- within(controlGLPart, isNeglectVPDEffect <- TRUE)
      resParms <- resParmsWithoutVPD <- partGLFitLRCWindows(dsR,
        nRecInDay = nRecInDay,
        dsTempSens = dsTempSens,
        controlGLPart = ctrlNeglectVPD,
        lrcFitter = lrcFitter
      )
      dsAns$FP_OPT_NoVPD[resParmsWithoutVPD$iCentralRec] <-
        resParmsWithoutVPD$resOpt
      dsAnsFluxes2 <- partGLInterpolateFluxes(dsR$Rg, 
        ds[[VPDVar]], ds[[TempVar]] # do prediction also using gap-Filled values
        , resParmsWithoutVPD,
        controlGLPart = ctrlNeglectVPD,
        lrcFitter = lrcFitter
      )
      dsAnsFluxes[iNAVPD, ] <- dsAnsFluxes2[iNAVPD, ]
    }
  }
  # Keenan et al., 2019
  if (isTRUE(controlGLPart$useNightimeBasalRespiration)) {
    dsAnsFluxes <- .computeRecoNight(dsAnsFluxes, dsR, dsTempSens)
  }

  nNAGPP <- sum(is.na(dsAnsFluxes$GPP))
  if (nNAGPP) warning("could not predict GPP in ", nNAGPP, " cases.")
  
  dsAns$FP_dRecPar <- dsAnsFluxes$dRecNextEstimate
  # copy quality flag from parameter row
  # 	first copy from iCentralRec to colNameAssoc
  # 	next transfer by dRecNextEstimate to each row
  iFiniteMeanRec <- which(is.finite(resParms[[colNameAssoc]]))
  dsAns$FP_qc[resParms[[colNameAssoc]][iFiniteMeanRec]] <-
    dsAns$FP_qc[resParms$iCentralRec[iFiniteMeanRec]]
  dsAns$FP_qc <- dsAns$FP_qc[1:nrow(dsAns) + dsAns$FP_dRecPar]
  # set quality flag to 2 where next parameter estimate is more than 14 days away
  dsAns$FP_qc[abs(dsAns$FP_dRecPar) > (14 * nRecInDay)] <- 2L
  # dsAns[is.finite(dsAns$FP_beta), ]
  
  units = attr(dsR$NEE, "units")
  RecoDTVar.s   <- paste0("Reco_DT", suffixDash)
  GPPDTVar.s    <- paste0("GPP_DT", suffixDash)

  dsAns[[RecoDTVar.s]] <- dsAnsFluxes$Reco %>% 
    set_attributes(list(varnames = RecoDTVar.s, units = units))
  
  dsAns[[GPPDTVar.s]] <- dsAnsFluxes$GPP %>% 
    set_attributes(list(varnames = GPPDTVar.s, units = units))
  
  if (controlGLPart$isSdPredComputed) {
    RecoDTSdVar.s <- paste0(RecoDTVar.s, "_SD")
    GPPDTSdVar.s <- paste0(GPPDTVar.s, "_SD")

    dsAns[[RecoDTSdVar.s]] <- dsAnsFluxes$sdReco %>% 
      set_attributes(list(varnames = RecoDTSdVar.s, units = units))

    dsAns[[GPPDTSdVar.s]] <- dsAnsFluxes$sdGPP %>% 
      set_attributes(list(varnames = GPPDTSdVar.s, units = units))
  }
  # sTEMP$GPP_DT_fqc <<- cbind(sDATA, sTEMP)[, QFFluxVar.s]
  # ! New code: MDS gap filling information are not copied from NEE_fmet
  # and NEE_fwin to GPP_fmet and GPP_fwin (since not known within this pure partitioning function)
  return(dsAns)
}

#' Relevant columns from original input with defined names
#' 
#' @param ds dataset with all the specified input columns and full days in equidistant times
#' @param NEEVar Variable of NEE
#' @param QFNEEVar Quality flag of variable
#' @param QFNEEValue Value of quality flag for _good_ (original) data
#' @param NEESdVar Variable of standard deviation of net ecosystem fluxes
#' 
#' @param TempVar Filled air or soil temperature variable (degC)
#' @param QFTempVar Quality flag of filled temperature variable
#' @param QFTempValue Value of temperature quality flag for _good_ (original) data
#' 
#' @param VPDVar Filled Vapor Pressure Deficit, VPD (hPa)
#' @param QFVPDVar Quality flag of filled VPD variable
#' @param QFVPDValue Value of VPD quality flag for _good_ (original) data
#' 
#' @param RadVar Filled radiation variable
#' @param QFRadVar Quality flag of filled radiation variable
#' @param QFRadValue Value of radiation quality flag for _good_ (original) data
#' 
#' @param PotRadVar Variable name of potential rad. (W / m2)
#' @param suffix string inserted into column names before identifier for NEE column defaults
#' 
#' @param controlGLPart further default parameters
#' 
#' @seealso [partGLControl()], [sEddyProc_sMDSGapFillAfterUstar()]
#' @export
partGLExtractStandardData <- function(ds
  , NEEVar = paste0("NEE", suffixDash, "_f")
  , QFNEEVar = if (!missing(QFNEEVar.s)) QFNEEVar.s else paste0("NEE", suffixDash, "_fqc")
  , QFNEEValue = if (!missing(QFNEEValue.n)) QFNEEValue.n else 0
  , NEESdVar = if (!missing(NEESdVar.s)) NEESdVar.s else paste0("NEE", suffixDash, "_fsd")
  , TempVar = paste0("Tair_f")
  , QFTempVar = if (!missing(QFTempVar.s)) QFTempVar.s else paste0("Tair_fqc")
  , QFTempValue = if (!missing(QFTempValue.n)) QFTempValue.n else 0
  , VPDVar = if (!missing(VPDVar.s)) VPDVar.s else paste0("VPD_f")
  , QFVPDVar = if (!missing(QFVPDVar.s)) QFVPDVar.s else paste0("VPD_fqc")
  , QFVPDValue = if (!missing(QFVPDValue.n)) QFVPDValue.n else 0
  , RadVar = if (!missing(RadVar.s)) RadVar.s else "Rg_f"
  , QFRadVar = if (!missing(QFRadVar.s)) QFRadVar.s else paste0("Rg_fqc")
  , QFRadValue = if (!missing(QFRadValue.n)) QFRadValue.n else 0
  , PotRadVar = if (!missing(PotRadVar.s)) PotRadVar.s else "PotRad_NEW"
  , suffix = if (!missing(Suffix.s)) Suffix.s else ""
  , NEEVar.s ## << deprecated
  , QFNEEVar.s ## << deprecated
  , QFNEEValue.n ## << deprecated
  , NEESdVar.s ## << deprecated
  , TempVar.s ## << deprecated
  , QFTempVar.s ## << deprecated
  , QFTempValue.n ## << deprecated
  , VPDVar.s ## << deprecated
  , QFVPDVar.s ## << deprecated
  , QFVPDValue.n ## << deprecated
  , RadVar.s ## << deprecated
  , QFRadVar.s ## << deprecated
  , QFRadValue.n ## << deprecated
  , PotRadVar.s ## << deprecated
  , Suffix.s ## << deprecated
  , controlGLPart = partGLControl() ## << further default parameters,
) {
  if (!missing(NEEVar.s)) NEEVar <- NEEVar.s # in default, lines too wide
  if (!missing(TempVar.s)) TempVar <- TempVar.s
  varNamesDepr <- c(
    "NEEVar.s", "TempVar.s", "VPDVar.s", "RadVar.s", "Suffix.s",
    "QFNEEVar.s", "QFNEEValue.n", "QFTempVar.s", "QFTempValue.n",
    "QFVPDVar.s", "QFVPDValue.n", "QFRadVar.s", "QFRadValue.n", "PotRadVar.s")
  varNamesNew <- c(
    "NEEVar", "TempVar", "VPDVar", "RadVar", "suffix",
    "QFNEEVar", "QFNEEValue", "QFTempVar", "QFTempValue",
    "QFVPDVar", "QFVPDValue", "QFRadVar", "QFRadValue", "PotRadVar")
  iDepr <- which(!c(
    missing(NEEVar.s), missing(TempVar.s), missing(VPDVar.s), missing(RadVar.s),
    missing(Suffix.s), missing(QFNEEVar.s), missing(QFNEEValue.n),
    missing(QFTempVar.s), missing(QFTempValue.n), missing(QFVPDVar.s),
    missing(QFVPDValue.n), missing(QFRadVar.s), missing(QFRadValue.n),
    missing(PotRadVar.s)
  ))
  if (length(iDepr)) {
    warning(
      "Argument names ", varNamesDepr[iDepr], " have been deprecated.",
      " Please, use instead ", varNamesNew[iDepr])
  }

  # Check if specified columns exist in sDATA or sTEMP and if numeric
  # and plausible. Then apply quality flag
  suffixDash <- paste0((if (fCheckValString(suffix)) "_" else ""), suffix)
  # Check if specified columns exist in sDATA or sTEMP and if numeric
  # and plausible. Then apply quality flag
  fCheckColNames(ds, c(
    NEEVar, QFNEEVar, TempVar, QFTempVar, RadVar, VPDVar,
    QFRadVar, PotRadVar, NEESdVar
  ), "sGLFluxPartition")
  fCheckColNum(ds, c(
    NEEVar, QFNEEVar, TempVar, QFTempVar, RadVar, VPDVar,
    QFRadVar, PotRadVar, NEESdVar
  ), "sGLFluxPartition")
  fCheckColPlausibility(ds, c(
    NEEVar, QFNEEVar, TempVar, QFTempVar, VPDVar,
    RadVar, QFRadVar, PotRadVar
  ), "sGLFluxPartition")
  NEEFiltered <- fSetQF(ds, NEEVar, QFNEEVar, QFNEEValue, "sGLFluxPartition")
  #
  # Apply quality flag for temperature and VPD
  # TODO: docu meteo filter, standard FALSE
  NEW_FP_Temp <- if (isTRUE(controlGLPart$isFilterMeteoQualityFlag)) {
    fSetQF(
      ds, TempVar, QFTempVar, QFTempValue, "partGLExtractStandardData")
  } else {
    ds[[TempVar]]
  }
  NEW_FP_VPD <- if (isTRUE(controlGLPart$isFilterMeteoQualityFlag)) {
    fSetQF(ds, VPDVar, QFVPDVar, QFVPDValue, "partGLExtractStandardData")
  } else {
    ds[[VPDVar]]
  }
  NEW_FP_Rg <- if (isTRUE(controlGLPart$isFilterMeteoQualityFlag)) {
    fSetQF(ds, RadVar, QFRadVar, QFRadValue, "partGLExtractStandardData")
  } else {
    ds[[RadVar]]
  }

  # Filter night time values only
  # ! Note: Rg <= 4 congruent with Lasslop et al., 2010 to define Night
  # for the calculation of E_0.n
  # Should be unfilled (original) radiation variable, therefore dataframe
  # set to sDATA only
  # ! New code: PotRad in sGLFluxPartition: Slightly different subset than
  # PV-Wave due to time zone correction (avoids timezone offset between Rg and PotRad)
  isPotRadZero <- if (isTRUE(controlGLPart$isNeglectPotRadForNight)) {
    TRUE
  } else {
    (ds[[PotRadVar]] <= 0)
  }
  isNotPotRadZero <- if (isTRUE(controlGLPart$isNeglectPotRadForNight)) {
    TRUE
  } else {
    (ds[[PotRadVar]] > 0)
  }
  isNight <- (ds[[RadVar]] <= 4 & isPotRadZero)
  # Filter day time values only
  # ! Note: Rg > 4 congruent with Lasslop et al., 2010 to define Day for the
  # calculation of paremeters of Light Response Curve
  # Should be unfilled (original) radiation variable, therefore dataframe set
  # to sDATA only, twutz does not understand this comment
  isDay <- (ds[[RadVar]] > 4 & isNotPotRadZero)
  #
  ## value<< a data.frame with columns
  dsR <- data.frame(
    sDateTime = ds[[1]] ## << first column of \code{ds},
    ## usually the time stamp
    ## not used, but usually first column is a DateTime is kept
    ## for aiding debug
    , NEE = NEEFiltered ## << NEE filtered for quality flay
    , sdNEE = ds[[NEESdVar]] ## << standard deviation of NEE
    ## with missing values replaced
    , Temp = NEW_FP_Temp ## << Temperature, quality filtered
    ## if isTRUE(controlGLPart$isFilterMeteoQualityFlag)
    , VPD = NEW_FP_VPD ## << Water pressure deficit, quality
    ## filtered if isTRUE(controlGLPart$isFilterMeteoQualityFlag)
    , Rg = NEW_FP_Rg ## << Incoming radiation
    , isDay = isDay ## << Flag that is true for daytime records
    , isNight = isNight ## << Flag that is true for nighttime records
  )
  ## details<<
  ## The LRC fit usually weights NEE records by its uncertainty. In order
  ## to also use
  ## records with missing \code{NEESdVar}, uncertainty of the missing values
  ## is by default set
  ## to a conservatively high value, parameterized by
  ## \code{controlGLPart$replaceMissingSdNEEParms)}.
  ## Controlled by argument \code{replaceMissingSdNEEParms} in
  ## \code{\link{partGLControl}}, but overruled
  ## by argument \code{neglectNEEUncertaintyOnMissing}.
  if (!controlGLPart$neglectNEEUncertaintyOnMissing) {
    dsR$sdNEE <- replaceMissingSdByPercentage(
      dsR$sdNEE, dsR$NEE,
      controlGLPart$replaceMissingSdNEEParms[1],
      controlGLPart$replaceMissingSdNEEParms[2]
    )
  }
  dsR
}


#' Estimate successive Rectangular Hyperbolic Light Response Curve parameters
#' 
#' @param ds data.frame with numeric columns NEE, sdNEE, Temp (degC), VPD, Rg,
#' and logical columns isNight and isDay
#' @param winSizeDays Window size in days for daytime fits
#' @param strideInDays step in days for shifting the windows
#' @param nRecInDay number of records within one day (for half-hourly data its 48)
#' @param dsTempSens data.frame that reports for each window temperature sensitivity
#' parameters E0 and RRef
#' @param isVerbose set to FALSE to suppress messages
#' @param controlGLPart list of further default parameters
#' @param lrcFitter R5 class instance responsible for fitting the light response curve
#' 
#' @return a tibble with a row for each window providing information on window size,
#' parameters estimated from light response curve (LRC), and their standard deviation.
#' Moreover, and estimated from night-time data after smoothing and forward-filling:
#' the uncertainty of temperature sensitivity \code{E0_night_sd}, And the respiration
#' at reference temperature \code{RRef_night}. Finally it contains a list column
#' \code{resOpt} with the optimization results.
#' 
#' @seealso [partGLFitLRCOneWindow()]
partGLFitLRCWindows <- function( 
  ds, winSizeDays = 4L, strideInDays = 2L, nRecInDay = 48L, dsTempSens, 
  isVerbose = TRUE, controlGLPart = partGLControl(), lrcFitter = RectangularLRCFitter()) {
  if (isVerbose) {
    message("  Estimating light response curve parameters from day time NEE ",
      appendLF = FALSE)
  }
  resLRC <- applyWindows(ds, partGLFitLRCOneWindow,
    prevRes = list(resOpt = list(thetaOpt = rep(NA_real_, 6L))),
    winSizeInDays = 4L,
    isVerbose = isVerbose,
    nRecInDay = nRecInDay
    
    , E0Win = dsTempSens,
    controlGLPart = controlGLPart,
    lrcFitter = lrcFitter
  )
  lrcSummary <- lapply(resLRC$resFUN, "[[", "summary")
  iNoSummary <- which(sapply(lrcSummary, length) == 0)
  if (length(iNoSummary)) {
    stop("expected summary returned by all fits, but found missing summaries.")
  }
  resOptList <- lapply(resLRC$resFUN, "[[", "resOpt")
  resOptDf <- tibble(resOpt = resOptList)
  resParms <- as_tibble(cbind(resLRC$winInfo, bind_rows(lrcSummary), resOptDf))

  # E0_night equals E0, but uncertainty might differ
  resParms$E0_bootstrap_sd <- resParms$E0_sd # due to bootstrap,
  # this may differ, save before overriding by night-time estimate
  resParms$E0 <- dsTempSens$E0 # report E0 and E0_sd from nighttime
  # even when LRC fit did not converge
  resParms$E0_sd <- dsTempSens$sdE0
  resParms$RRef_night <- dsTempSens$RRef
  # summary$iMeanRec yet based on window instead of entire time, need to
  # add beginning of window
  resParms$iMeanRec <- resParms$iRecStart - 1L + resParms$iMeanRec
  # 	# omit records where NULL was returned
  # 	iWinNoFit <- which(is.na(resParms$parms_out_range) )
  # 	if (length(iWinNoFit) ) {
  # 		resParms <- resParms[-iWinNoFit, ]
  # 		resParms$resOptList <- resParms$resOptList[-iWinNoFit]
  # 	}
  ## value<< a tibble with a row for each window
  ## providing information on window size, parameters estimated
  ## from light response curve (LRC),
  ## and their standard deviation.
  ## Moreover, and estimated from night-time data after smoothing
  ## and forward-filling:
  ## the uncertainty of temperature sensitivity \code{E0_night_sd},
  ## And the respiration at reference temperature \code{RRef_night}.
  ## Finally it contains a list column\code{resOpt} with the optimization results.
  resParms
}

#' Estimate Rectangular Hyperbolic Light Response Curve parameters
#' 
#' @details Estimation as in Lasslop et al., 2010 for successive periods, i.e.
#' windows. Fitting is done on a subset of the data where isDay, NEE, sdNEE,
#' Temp, Rg, and VPD are all non-NA and isDay is TRUE
#' 
#' @param ds data.frame with numeric columns NEE, sdNEE, Temp (degC), VPD, Rg,
#' and logical columns isNight and isDay
#' @param winInfo one-row data.frame with window information, including iWindow
#' @param prevRes component prevRes from previous result, here with item prevE0
#' @param E0Win data.frame with columns E0, sdE0, RRef from nighttime,
#' one row for each window
#' @param controlGLPart list of further default parameters
#' @param lrcFitter R5 class instance responsible for fitting the light response curve
#' 
#' @return a list with components resOpt, summary, isValid
#' @seealso [partitionNEEGL()]
partGLFitLRCOneWindow <- function(
  ds, winInfo, prevRes, E0Win, controlGLPart = partGLControl(), lrcFitter = RectangularLRCFitter()) {
  
  isValidDayRecNoSdNEEConstraint <- !is.na(ds$isDay) & ds$isDay &
    !is.na(ds$NEE) & !is.na(ds$Temp) & !is.na(ds$Rg)
  ## details<<
  ## With option \code{controlGLPart$neglectNEEUncertaintyOnMissing = TRUE}
  ## all sdNEE are set to 1 if there is any missing sdNEE at
  ## otherwise valid records
  if (isTRUE(controlGLPart$neglectNEEUncertaintyOnMissing) &&
    any(!is.finite(ds$sdNEE[isValidDayRecNoSdNEEConstraint]))) {
    ds$sdNEE <- 1L
  }
  isValidDayRecNoVPDConstraint <- isValidDayRecNoSdNEEConstraint &
    !is.na(ds$sdNEE)
  ## details<<
  ## If  \code{controlGLPart$isNeglectVPDEffect = TRUE}, also records
  ## with VPD = NA maybe valid
  isValidDayRec <- if (isTRUE(controlGLPart$isNeglectVPDEffect)) {
    isValidDayRecNoVPDConstraint
  } else {
    isValidDayRecNoVPDConstraint & !is.na(ds$VPD)
  }
  iMeanRecInDayWindow <- as.integer(round(mean(which(isValidDayRec))))
  
  # if too few records or if no temperature-respiration relationship could be
  # found, indicate no-fit, but report Window properties
  E0 <- E0Win$E0[winInfo$iWindow]
  # when provided fixed E0 the column sdEO does not exist
  sdE0 <- if (length(E0Win$sdE0)) E0Win$sdE0[winInfo$iWindow] else NA
  RRefNight <- E0Win$RRef[winInfo$iWindow]

  getNAResult <- function(convergenceCode) {
    list(
      resOpt = NULL,
      summary = data.frame(
        nValidRec = sum(isValidDayRec),
        iMeanRec = iMeanRecInDayWindow,
        convergence = convergenceCode,
        E0 = E0, E0_sd = sdE0,
        RRefNight = RRefNight
      ),
      isValid = FALSE)
  }
  if (is.na(E0)) return(getNAResult(1010L))
  
  # details<< If there are too few records (<
  # \code{controlGLPart$minNRecInDayWindow}), then the constraint on non-NA VPD
  # is neglected and \code{controlGLPart$isNeglectVPDEffect} is set to TRUE If
  # there are still too few records, an NA-result with convergence code 1011L is
  # returned.
  if ((sum(isValidDayRec) < controlGLPart$minNRecInDayWindow)) {
    controlGLPart$isNeglectVPDEffect <- TRUE
    isValidDayRec <- isValidDayRecNoVPDConstraint
    iMeanRecInDayWindow <- as.integer(round(mean(which(isValidDayRec))))
    if ((sum(isValidDayRec) < controlGLPart$minNRecInDayWindow)) {
      return(getNAResult(1011L))
    }
  }
  dsDay <- ds[isValidDayRec, ]
  ## details<<
  ## Each window estimate is associated with a time or equivalently with a record.
  ## The first record, i.e. row number, of the day-window is reported.
  ## Moreover, the mean of all valid records numbers in the daytime window is
  ## reported for interpolation.
  # better report NA and care for it properly: if (is.na(iMeanRecInDayWindow))
  # iMeanRecInDayWindow <- as.integer(nrow(ds)%/%2)
  # TODO firstRecInDayWindow.i <- which(SubsetDayPeriod.b)[1]
  # the rownumber of the first record inside the day window
  ## seealso<< \code{\link{partGL_FitNight_1win_E0_kernel}}
  #
  ## seealso<< \code{\link{LRC_fitLRC}}
  # lrcFitter <- RectangularLRCFitter()
  resOpt <- resOpt0 <- lrcFitter$fitLRC(dsDay,
    E0 = E0, sdE0 = sdE0,
    RRefNight = RRefNight,
    controlGLPart = controlGLPart,
    lastGoodParameters = prevRes$resOpt$thetaOpt
  )
  if (!is.finite(resOpt$thetaOpt[1])) return(getNAResult(resOpt$convergence))
  
  sdTheta <- resOpt$thetaOpt
  sdTheta[] <- NA
  sdTheta[resOpt$iOpt] <- sqrt(diag(resOpt$covParms)[resOpt$iOpt])
  
  # compute GPP2000 without VPD effect
  GPP2000 <- lrcFitter$predictLRC(
    resOpt$thetaOpt, 2000, VPD = NA, Temp = NA, fixVPD = TRUE)$GPP
  ans <- list(
    resOpt = resOpt,
    summary = cbind(
      data.frame(
        nValidRec = nrow(dsDay),
        iMeanRec = iMeanRecInDayWindow,
        convergence = resOpt$convergence,
        parms_out_range = as.integer(!all(1:5 %in% resOpt$iOpt))
      ),
      as.data.frame(t(resOpt$thetaOpt)),
      as.data.frame(t(structure(sdTheta, names = paste0(names(sdTheta), "_sd")))),
      GPP2000 = GPP2000, isValid = TRUE
    )
  )
  return(ans)
}

.bootStrapLRCFit <- function(
  ### Compute parameters uncertainty by bootstrap
  theta0, iOpt, dsDay, sdE_0.n, parameterPrior, controlGLPart,
  lrcFitter ## << Light Response Curve R5 instance
  , iPosE0 = 5L ## << position (integer scalar) of temperature
  ## sensitivity in parameter vector
  # should be dealt with iOpt: isNeglectVPD = FALSE
) {
  ## value<<
  ## Matrix with each row a parameter estimate on a different bootstrap sample.
  ## Some of the rows might be NA, if the subsampled dataset could not be fitted.
  ans <- matrix(NA,
    nrow = controlGLPart$nBootUncertainty, ncol = length(theta0),
    dimnames = list(NULL, names(theta0))
  )
  ## details<<
  ## In addition to resampling the original data, also the temperature sensitivity
  ## is resampled
  ## from its uncertainty distribution but bounded in between 50 to 400.
  E0r <- rnorm(controlGLPart$nBootUncertainty, theta0[5L], sdE_0.n)
  E0 <- pmax(50, pmin(400, E0r))
  theta <- theta0
  # iBoot <- 1L
  for (iBoot in c(1:controlGLPart$nBootUncertainty)) {
    idx <- sample(nrow(dsDay), replace = TRUE)
    dsDayB <- dsDay[idx, ]
    theta[iPosE0] <- E0[iBoot]
    resOptBoot <- try(lrcFitter$optimLRCOnAdjustedPrior(theta,
      iOpt = iOpt,
      dsDay = dsDayB,
      parameterPrior = parameterPrior, ctrl = controlGLPart
    ), silent = TRUE)
    if (!inherits(resOptBoot, "try-error") && resOptBoot$convergence == 0L) {
      # TODO: also remove the very bad cases?
      ans[iBoot, ] <- resOptBoot$theta
    } else {
      # recover()
    }
  }
  ans
}

#' Interpolate Reco and GPP from two neighboring parameter sets
#' 
#' @details `resLRC$iFirstRecInCentralDay` must denote the row for which the
#' LRC parameters are representative, here, the first record of the center day
#' 
#' 
#' @param Rg numeric vector of photosynthetic flux density [umol / m2 / s] or Global Radiation
#' @param VPD numeric vector of Vapor Pressure Deficit [hPa]
#' @param Temp numeric vector of Temperature [degC]
#' @param resParms data frame with results of \code{partGLFitLRCWindows} of
#' fitting the light-response-curve for several windows
#' @param controlGLPart further default parameters, see \code{\link{partGLControl}}
#' @param lrcFitter R5 class instance responsible for fitting the LRC curve
#' @param isVerbose set to FALSE to suppress messages
#' 
#' @seealso [partitionNEEGL()], [partGLFitLRCWindows()]
partGLInterpolateFluxes <- function(
  Rg, VPD, Temp, resParms, controlGLPart = partGLControl(), 
  lrcFitter, isVerbose = TRUE) {
  # create a dataframe with index of rows of estimates before and after and
  # corresponding weights
  iValidWin <- which(is.finite(resParms$parms_out_range))
  summaryLRC <- resParms %>% select(!"resOpt") %>% slice(iValidWin)
  
  resOptList <- resParms$resOpt[iValidWin]
  nLRC <- nrow(summaryLRC)
  nRec <- length(Rg)
  Temp_Kelvin <- Temp + 273.15
  if (isTRUE(controlGLPart$isAssociateParmsToMeanOfValids)) {
    # there might be several rows with the same iMeanRec, omit those rows
    # unless the first of each reoccuring iMeanRec
    # they all will be based on the same valid data, if iMeanRec is equal
    tabMeanRec <- table(summaryLRC$iMeanRec)
    iRecsDouble <- as.integer(names(tabMeanRec[tabMeanRec > 1L]))
    iRecsOmit <- do.call(c, lapply(iRecsDouble, function(iRecDouble) {
      which(summaryLRC$iMeanRec == iRecDouble)[-1L]
    }))
    if (length(iRecsOmit)) summaryLRC <- summaryLRC[-iRecsOmit, ]
  }
  ## details<<
  ## Parameter estimates are reported for the central record of the window, or
  ## if \code{isTRUE(controlGLPart.l$isAssociateParmsToMeanOfValids)} for the
  ## mean time of all valid records within the window
  # for each original record merge parameters assicated with previous fit or
  # next fit respectively
  colNameAssoc <- ifelse(isTRUE(controlGLPart$isAssociateParmsToMeanOfValids), 
    "iMeanRec", "iCentralRec")
  dsAssoc <- .partGPAssociateSpecialRows(summaryLRC[[colNameAssoc]], nRec)
  # now we have columns iBefore and iAfter, which can be used to merge the
  # parameter estimate columns to each row
  parNames <- lrcFitter$getParamNames()
  dsBefore <- merge(
    structure(data.frame(dsAssoc$iSpecialBefore, dsAssoc$iBefore),
      names = c("iParRec", colNameAssoc)),
    summaryLRC[, c(colNameAssoc, parNames)]
  )
  dsAfter <- merge(
    structure(data.frame(dsAssoc$iSpecialAfter, dsAssoc$iAfter),
      names = c("iParRec", colNameAssoc)),
    summaryLRC[, c(colNameAssoc, lrcFitter$getParamNames())]
  )
  if ((nrow(dsBefore) != nRec) || (nrow(dsAfter) != nRec)) {
    stop("error in merging parameters to original records.")
  }
  Reco2 <- lapply(list(dsBefore, dsAfter), function(dsi) {
    # twutz170316: fLloydTaylor gives unreasonable values with very
    # low temperatures, hence constrain lower temperature
    tmp <- fLloydTaylor(dsi$RRef, dsi$E0, pmax(-40, Temp_Kelvin),
      TRef = 273.15 + 15)
  })
  # dsi <- dsBefore
  GPP2 <- lapply(list(dsBefore, dsAfter), function(dsi) {
    theta <- as.matrix(dsi[, parNames])
    tmp <- lrcFitter$predictLRC(theta, Rg, VPD, Temp = pmax(-40, Temp))$GPP
  })
  # interpolate between previous and next fit, weights already sum to 1
  Reco <- (dsAssoc$wBefore * Reco2[[1]] + dsAssoc$wAfter * Reco2[[2]])
  # / (dsAssoc$wBefore + dsAssoc$wAfter)
  GPP <- (dsAssoc$wBefore * GPP2[[1]] + dsAssoc$wAfter * GPP2[[2]])
  # / (dsAssoc$wBefore + dsAssoc$wAfter)
  ans <- ansPred <- data.frame(Reco, GPP)
  if (isTRUE(controlGLPart$isSdPredComputed)) {
    # dsi <- dsBefore
    varPred2 <- lapply(list(dsBefore, dsAfter), function(dsi) {
      theta <- as.matrix(dsi[, parNames])
      grad <- lrcFitter$computeLRCGradient(theta, Rg, VPD, pmax(-40, Temp))
      varPred <- matrix(NA_real_,
        nrow = nrow(dsi), ncol = 2L,
        dimnames = list(NULL, c("varGPP", "varReco")))
      # iRec <- 1L
      for (iRec in 1:nrow(dsi)) {
        iParRec <- dsi$iParRec[iRec] # row within sequence of valid
        # parameters (at centralREc or meanRec)
        # get the fitting object, TODO better document
        resOpt <- resOptList[[iParRec]]
        gradGPP <- grad$GPP[iRec, ]
        gradReco <- grad$Reco[iRec, ]
        # make sure parameter names match postions in covParms
        varPred[iRec, 2L] <- varReco <- gradReco %*%
          resOpt$covParms[names(gradReco), names(gradReco)] %*% gradReco
        varPred[iRec, 1L] <- varGPP <- gradGPP %*%
          resOpt$covParms[names(gradGPP), names(gradGPP)] %*% gradGPP
      }
      varPred
    })
    sdGPP <- sqrt(dsAssoc$wBefore^2 * varPred2[[1]][, "varGPP"] +
      dsAssoc$wAfter^2 * varPred2[[1]][, "varGPP"])
    sdReco <- sqrt(dsAssoc$wBefore^2 * varPred2[[1]][, "varReco"] +
      dsAssoc$wAfter^2 * varPred2[[1]][, "varReco"])
    ans <- cbind(ansPred, data.frame(sdGPP = sdGPP, sdReco = sdReco))
  }
  # compute differences in rows to next parameter estimate
  dRecBefore <- dsBefore[[colNameAssoc]] - 1:nrow(dsBefore)
  dRecAfter <- dsAfter[[colNameAssoc]] - 1:nrow(dsBefore)
  isBeforeCloser <- abs(dRecBefore) <= abs(dRecAfter)
  ans$dRecNextEstimate <- ifelse(isBeforeCloser, dRecBefore, dRecAfter)
  ## value<< data.frame with nrow() rows and columns  GPP, Reco, varGPP
  ## , varReco, and dRecNextEstimate
  ans
}

#' associate each row with the previous and next row from a subset of rows
#' @param iRowsSpecial ordered unique integer vector specifying the rows for which
#'  some special data is available
#' @param nRec integer scalar of number of rows in the full data.frame
#' 
#' @return a dataframe with index of previous and next rows inside the subset
.partGPAssociateSpecialRows <- function(iRowsSpecial, nRec) {
  ## details<<
  ## When only for a subset of rows some more data available, i.e parameter
  ## estimates for
  ## only a subset of rows, this function creates
  ## columns that refer to the previous and next row that are in the subset.
  ## E.g. if some more data is available for rows 3 and 7, then rows 4:6 will
  ## indicate\code{iBefore = 3, iAfter = 7}.
  ans <- data.frame(
    iRec = 1:nRec ## << the original row number
    , iSpecialBefore = NA_integer_ ## << index within \code{iRowsSpecial}
    , iSpecialAfter = NA_integer_ ## << index within \code{iRowsSpecial}
    , iBefore = NA_integer_ ## << index of the previous special row
    , iAfter = NA_integer_ ## << index of the next special row
    , wBefore = NA_real_ ## << weight of the previous, inverse of the
    ## distance in records
    , wAfter = NA_real_
  ) ## << weight of the next, inverse of the
  ## distance in records
  ## details<<
  ## The subset rows inside the subset refer both (before and after) to the
  ## same subset rows, with weights 1 / 2
  nRecS <- length(iRowsSpecial)
  if (0 == nRecS) {
    stop("cannot associate special rows, if length of argument iRowsSpecial is zero.")
  }
  ans[iRowsSpecial, "iSpecialBefore"] <- ans[iRowsSpecial, "iSpecialAfter"] <-
    seq_along(iRowsSpecial)
  ans[iRowsSpecial, "iBefore"] <- ans[iRowsSpecial, "iAfter"] <- iRowsSpecial
  ans[iRowsSpecial, c("wBefore", "wAfter")] <- 0.5
  # iS <- 2L
  for (iS in 1:nRecS) {
    currRec <- iRowsSpecial[iS]
    # before and after last special row will be treated afterwards
    prevRec <- if (iS == 1L) currRec else iRowsSpecial[iS - 1L]
    nextRec <- if (iS == nRecS) currRec else iRowsSpecial[iS + 1L]
    # c(prevRec, currRec, nextRec)
    ## details<<
    ## The weight is inversely proportional to the distance in rows
    ## The two weights wBefore and wAfter always sum to 1
    distPrev <- currRec - prevRec
    if (distPrev > 1L) {
      ans[(prevRec + 1L):(currRec - 1L), "iSpecialAfter"] <- iS
      ans[(prevRec + 1L):(currRec - 1L), "iAfter"] <- currRec
      ans[(prevRec + 1L):(currRec - 1L), "wAfter"] <- (1:(distPrev - 1)) / distPrev
    }
    distNext <- nextRec - currRec
    if (distNext > 1L) {
      ans[(currRec + 1L):(nextRec - 1L), "iSpecialBefore"] <- iS
      ans[(currRec + 1L):(nextRec - 1L), "iBefore"] <- currRec
      ans[(currRec + 1L):(nextRec - 1L), "wBefore"] <- ((distNext - 1):1) / distNext
    }
  }
  ## details<<
  ## the rows before the first subset row refer bot (after and before) to the
  ## first subset row with weights 1 / 2
  ## similar the rows after the last subset row refer to the last subset row
  ans[1:iRowsSpecial[1], c("iSpecialBefore", "iSpecialAfter")] <- 1L
  ans[1:iRowsSpecial[1], c("iBefore", "iAfter")] <- iRowsSpecial[1L]
  ans[1:iRowsSpecial[1], c("wBefore", "wAfter")] <- 0.5
  ans[iRowsSpecial[nRecS]:nrow(ans), c("iSpecialBefore", "iSpecialAfter")] <- nRecS
  ans[iRowsSpecial[nRecS]:nrow(ans), c("iBefore", "iAfter")] <- iRowsSpecial[nRecS]
  ans[iRowsSpecial[nRecS]:nrow(ans), c("wBefore", "wAfter")] <- 0.5
  ans
}


#' recompute ecosystem respiration at nighttime with nighttime reference
#' 
#' @details If `partGLControl$useNightimeBasalRespiration == TRUE` then
#' nighttime respiration is recomputed based on basal respiration inferred
#' from nighttime data.
#' 
#' @param dsAnsFluxes data.frame with columns Reco and sdReco
#' @param dsR data.frame with columns isNight, and Temp
#' @param dsTempSens data.frame with column `iCentralRec`, `E0`, and `RRef`
#' 
#' @return data.frame with updated night-time `Reco` and `sdReco`
.computeRecoNight <- function(dsAnsFluxes, dsR, dsTempSens) {
  # Calculate the ecosystem respiration Reco by LlyodAndTaler
  TRefK <- 273.15 + 15
  T0 <- 227.13
  RRef <- approx(dsTempSens$iCentralRec, dsTempSens$RRef, 1:nrow(dsR), rule = 2)$y
  E0 <- approx(dsTempSens$iCentralRec, dsTempSens$E0, 1:nrow(dsR), rule = 2)$y

  inds = dsR$isNight
  dsAnsFluxes$Reco[inds] <- respNight <- 
    fLloydTaylor(RRef[inds], E0[inds], fConvertCtoK(dsR$Temp[inds]), TRef = TRefK)
  ## details<< With given uncertainty of temperature sensitivity E0, the
  ## logarithm of respraiton R is normally distributed, i.e. R is lognormally
  ## distributed. With low uncertainty (sigma*) < 1.2 this is well
  ## approximated by a normal distribution. Hence the sqrt(second momemt)
  ## of the lognormal is reported as sdReco
  # log(R) = log(RRef) + tempFac*E0
  sdE0 <- sqrt(approx(dsTempSens$iCentralRec, dsTempSens$sdE0^2, 1:nrow(dsR), rule = 2)$y)
  tempFacLloydTaylor <- 1 / (TRefK - T0) - 1 / (fConvertCtoK(dsR$Temp) - T0)
  sigmaLogR <- (abs(tempFacLloydTaylor) * sdE0)[inds]
  
  # formula of variance of lognormal with mu replaced by mean
  # see lognrom/inst/docu/varianceBySigmaAndExpected.Rmd
  dsAnsFluxes$sdReco[inds] <- sqrt(exp(sigmaLogR^2) - 1) * respNight
  dsAnsFluxes
}

# .tmp.f <- function() { # omit declaration for now
#   computeAggregatedCovariance <- function(
#     ### Covariances between Reco and GPP predictions due same uncertain model coefficients
#     dsPred ## << data.frame with predictors (Rg, VPD, Temp)
#     , resParms ## << data.frame with results of \code{partGLFitLRCWindows} of
#     ## fitting the light-response-curve for several windows
#     , resParmsNoVPD ## << data.frame like resParms, but was fitted with option
#     ## isNeglectVPDEffect = TRUE, for predicting if VPD is missing
#     , controlGLPart = partGLControl() ## << further default parameters,
#     ## see \code{\link{partGLControl}} with entry "isAssociateParmsToMeanOfValids"
#     , lrcFitter ## << R5 class instance responsible for fitting the light response
#     ## curve, with method getParamNames()
#     , iAggregate = 1:nrow(dsPred) ## << row indices about which to
#     ## sum over, must be contiguous
#   ) {
#     sumCovGPP <- 0
#     sumCovReco <- 0
#     # merge parameters to predictors
#     isValidWin <- is.finite(resParms$parms_out_range)
#     summaryLRC <- resParms[isValidWin, , drop = FALSE]
#     resOptList <- resParms$resOptList[isValidWin]
#     colNameAssoc <- if (isTRUE(controlGLPart$isAssociateParmsToMeanOfValids)) {
#       "iMeanRec"
#     } else {
#       "iCentralRec"
#     }
#     dsAssoc <- .partGPAssociateSpecialRows(summaryLRC[[colNameAssoc]], nRec)
#     # now we have columns iBefore and iAfter, which can be used to
#     # merge the parameter estimate columns to each row
#     parNames <- lrcFitter$getParamNames()
#     dsBefore <- merge(
#       structure(data.frame(dsAssoc$iSpecialBefore, dsAssoc$iBefore),
#         names = c("iParRec", colNameAssoc)),
#       summaryLRC[, c(colNameAssoc, parNames)]
#     )
#     dsAfter <- merge(
#       structure(data.frame(dsAssoc$iSpecialAfter, dsAssoc$iAfter),
#         names = c("iParRec", colNameAssoc)),
#       summaryLRC[, c(colNameAssoc, lrcFitter$getParamNames())]
#     )
#     #
#     # compute gradients (each entry is a list with components vectors Reco and GPP
#     grad2 <- lapply(list(dsBefore, dsAfter), function(dsi) {
#       theta <- as.matrix(dsi[, parNames])
#       grad <- lrcFitter$computeLRCGradient(
#         theta, dsPred$Rg, dsPred$VPD, pmax(-40, dsPred$Temp))
#     })
#     for (iRec in iAggregate[-1]) { # in first line only entry with i == j
#       iParRecBefore <- dsAssoc$iSpecialBefore[iRec]
#       resOptBefore <- resOptList[[iParRecBefore]]
#       iParRecAfter <- dsAssoc$iSpecialAfter[iRec]
#       resOptAfter <- resOptList[[iParRecAfter]]
#       gradGPPIBefore <- grad[[1]]$GPP[iRec, ]
#       gradRecoIBefore <- grad[[1]]$Reco[iRec, ]
#       gradGPPIAfter <- grad[[2]]$GPP[iRec, ]
#       gradRecoIAfter <- grad[[2]]$Reco[iRec, ]
#       for (jRec in iAggregate[1]:(iRec - 1)) { # iter from start to one before iRec
#         # assume covariance larger than zero only in same interval between estimates
#         if ((dsAssoc$iSpecialBefore[jRec] == iParRecBefore) &&
#           (dsAssoc$iSpecialAfter[jRec] == iParRecAfter) ) {
#           gradGPPJBefore <- grad[[1]]$GPP[jRec, ]
#           gradRecoJBefore <- grad[[1]]$Reco[jRec, ]
#           gradGPPJAfter <- grad[[2]]$GPP[jRec, ]
#           gradRecoJAfter <- grad[[2]]$Reco[jRec, ]
#           covRecoBefore <- gradRecoIBefore %*%
#             resOptBefore$covParms[names(gradRecoIBefore), names(gradRecoIBefore)] %*%
#             gradRecoJBefore
#           covGPPBefore <- gradGPPIBefore %*%
#             resOptBefore$covParms[names(gradGPPIBefore), names(gradGPPIBefore)] %*%
#             gradGPPJBefore
#           covRecoAfter <- gradRecoIAfter %*%
#             resOptAfter$covParms[names(gradRecoIAfter), names(gradRecoIAfter)] %*%
#             gradRecoJAfter
#           covGPPBAfter <- gradGPPIAfter %*%
#             resOptAfter$covParms[names(gradGPPIAfter), names(gradGPPIAfter)] %*%
#             gradGPPJAfter
#           covRecoMix <- dsAssoc$wBefore[iRec] * dsAssoc$wBefore[jRec] * covRecoBefore +
#             dsAssoc$wAfter[iRec] * dsAssoc$wAfter[jRec] * covRecoAfter
#           covGPPMix <- dsAssoc$wBefore[iRec] * dsAssoc$wBefore[jRec] * covGPPBefore +
#             dsAssoc$wAfter[iRec] * dsAssoc$wAfter[jRec] * covGPPAfter
#           sumCovReco <- sumCovReco + 2 * covRecoMix # symm, add also from transposed
#           sumCovGPP <- sumCovGPP + 2 * covGPPMix
#         } # match of iParRecBefore and iParRecAfter
#       } # for jRec
#     } # for iRec
#     ## value<< named numeri vector with two components:
#     ans <- c(
#       sumCovGPP = sumCovGPP ## << sum_{i<>j} cov(GPP_i, GPP_j)
#       , sumCovReco = sumCovReco ## << sum_{i<>j} cov(Reco_i, Reco_j)
#     )
#   }
# } # .tmp.f

# .tmp.f <- function() {
#   # execute on recover jsut before ans <- in .partGPInterploateFluxes,
#   # from test_that("interpolate Fluxes"
#   # there is only one period with changing parametrization
#   x <- seq_along(Rg)
#   plot(Reco ~ x)
#   lines(Reco2[[1]] ~ x, col = "blue")
#   lines(Reco2[[2]] ~ x, col = "orange")
#   # note how the black dots first follow the blue line, then converge
#   # to the orange line
#   plot(GPP ~ x)
#   lines(GPP2[[1]] ~ x, col = "blue")
#   lines(GPP2[[2]] ~ x, col = "orange")
#   # not much difference between blue and orange
#   # (difference only seen if plot window is very large)
# }
