## Unit tests for fConvertTimeToPosix functions +++
# Author: TW
context("LightResponseCurveFitter")

# 8 first days of June from IT-MBo.2005.txt
# 10 days from June from Example_DETha98.txt shipped with REddyProc
.tmp.f <- function() {
  # save(dsNEE, file="tmp/dsNEE_Tharandt.RData")
  load("tmp/dsNEE_Tharandt.RData") # dsNEE
  dsNEE$Temp <- dsNEE$Tair_f
  dsNEE$Rg_f <- dsNEE$Rg
  dsNEE$isNight <- (dsNEE$Rg_f <= 4 & dsNEE$PotRad_NEW == 0)
  dsNEE$isDay <- (dsNEE$Rg_f > 4 & dsNEE$PotRad_NEW != 0)
}

.tmp.f2 <- function() {
  # else stop in partitionNEEGL, and grap ds:
  attr(ds$sDateTime, "tzone") <- "UTC"
  dsJune <- ds
  dsNEE <- dsJune[
    dsJune$sDateTime >= as.POSIXct("1998-06-01", tz = "UTC") & dsJune$sDateTime < as.POSIXct("1998-06-10", tz = "UTC"),
    c("sDateTime", "NEE_f", "NEE_fqc", "NEE_fsd", "Tair_f", "Tair_fqc", "VPD_f", "VPD_fqc", "Rg_f", "PotRad_NEW")
  ]
  save(dsNEE, file = "tmp/dsNEE_Tharandt.RData")
  dput(dsNEE)
}

# test_file("tests/testthat/test_LightResponseCurveFitter.R")
source("data-LRC.R")

RectangularLRCFitter$methods()
lrcFitter <- RectangularLRCFitter()

test_predictLRC <- function(lrcFitter) {
  theta0 <- structure(c(0, 27.3333395589509, 0.162207578338878, 2.59392002410639, 185, 1.1), 
    .Names = c("k", "beta", "alpha", "RRef", "E0", "logitconv"))
  nRow <- 5L
  iRowK <- 4L
  theta0M <- do.call(rbind, lapply(1:nRow, \(i) theta0 ))
  theta0M[iRowK, "k"] <- 0.05
  dss <- subset(dsNEE, mday(sDateTime) %in% 1:8)
  dssDay <- subset(dss, isDay == TRUE)[rep(1, nRow), ]
  dsDay <- rename(dssDay, NEE = NEE_f, sdNEE = NEE_fsd, Rg = Rg_f, 
    Temp = Temp, VPD = VPD_f)
  
  # lrcFitter$trace("predictLRC", browser)   # lrcFitter$untrace("predictLRC")
  gpp <- gpp0 <- lrcFitter$predictLRC(theta0M, Rg = dsDay$Rg, VPD = dsDay$VPD, Temp = dsDay$Temp, VPD0 = 0)$GPP
  expect_equal(length(gpp), nRow)
  expect_true(all(gpp[-iRowK] == gpp[1]))
  expect_true(all(gpp[iRowK] != gpp[1]))
  #
  gpp <- lrcFitter$predictLRC(theta0M, Rg = dsDay$Rg, VPD = dsDay$VPD, Temp = dsDay$Temp, VPD0 = 0, fixVPD = FALSE)$GPP
  expect_equal(length(gpp), nRow)
  expect_true(all(gpp[-iRowK] == gpp[1])) # the only row with k!=0
  expect_true(all(gpp[iRowK] != gpp[1]))
  #
  gpp <- lrcFitter$predictLRC(theta0M, Rg = dsDay$Rg, VPD = rep(NA_real_, nRow), Temp = dsDay$Temp, VPD0 = 0, fixVPD = FALSE)$GPP
  expect_equal(length(gpp), nRow)
  expect_true(all(is.na(gpp))) # VPD in computation
  #
  gpp <- lrcFitter$predictLRC(theta0M, Rg = dsDay$Rg, VPD = rep(NA_real_, nRow), Temp = dsDay$Temp, VPD0 = 0, fixVPD = TRUE)$GPP
  expect_equal(length(gpp), nRow)
  expect_true(all(gpp == gpp[1])) # same value not regarding different k
  expect_true(is.finite(gpp[1])) # finite despite NA VPD
  #
  gpp <- lrcFitter$predictLRC(theta0M, Rg = dsDay$Rg, VPD = rep(NA_real_, nRow), Temp = dsDay$Temp, VPD0 = 0)$GPP
  expect_equal(length(gpp), nRow)
  expect_true(all(gpp[-iRowK] == gpp[1])) # the only row with k!=0
  expect_true(is.na(gpp[iRowK]))
}

test_that("RectangularLRCFitter$predictLRC vector form", {
  lrcFitter <- RectangularLRCFitter()
  test_predictLRC(lrcFitter)
  # lrcFitter <- LogisticSigmoidLRCFitter()
})

test_that("NonrectangularLRCFitter$predictLRC vector form", {
  lrcFitter <- NonrectangularLRCFitter()
  test_predictLRC(lrcFitter)
})

test_LRCGRadient_vectorForm <- function(lrcFitter) {
  theta0 <- structure(c(0, 27.3333395589509, 0.162207578338878, 2.59392002410639, 185, 1.1), .Names = c("k", "beta", "alpha", "RRef", "E0", "logitconv"))
  nRow <- 5L
  iRowK <- 4L
  theta0M <- do.call(rbind, lapply(1:nRow, function(i) {
    theta0
  }))
  theta0M[iRowK, "k"] <- 0.05
  dss <- subset(dsNEE, as.POSIXlt(dsNEE$sDateTime)$mday %in% 1:8)
  dssDay <- do.call(rbind, lapply(1:nRow, function(i) {
    subset(dss, isDay == TRUE)[1, ]
  }))
  dsDay <- data.frame(NEE = dssDay$NEE_f, sdNEE = dssDay$NEE_fsd, Rg = dssDay$Rg_f, Temp = dssDay$Temp, VPD = dssDay$VPD_f)
  #
  gpp <- gpp0 <- lrcFitter$computeLRCGradient(theta0M, Rg = dsDay$Rg, VPD = dsDay$VPD, Temp = dsDay$Temp, VPD0 = 0)$GPP
  expect_equal(nrow(gpp), nRow)
  expect_true(all(gpp[-iRowK, ][2, ] == gpp[1, ]))
  expect_true(all(gpp[iRowK, "k"] != gpp[1, "k"]))
  #
  gpp <- lrcFitter$computeLRCGradient(theta0M, Rg = dsDay$Rg, VPD = dsDay$VPD, Temp = dsDay$Temp, VPD0 = 0, fixVPD = FALSE)$GPP
  expect_equal(nrow(gpp), nRow)
  expect_true(all(gpp[-iRowK, ][2, ] == gpp[1, ]))
  expect_true(all(gpp[iRowK, "k"] != gpp[1, "k"]))
  #
  gpp <- lrcFitter$computeLRCGradient(theta0M, Rg = dsDay$Rg, VPD = rep(NA_real_, nRow), Temp = dsDay$Temp, VPD0 = 0, fixVPD = FALSE)$GPP
  expect_equal(nrow(gpp), nRow)
  expect_true(all(is.na(gpp))) # VPD in computation
  #
  gpp <- lrcFitter$computeLRCGradient(theta0M, Rg = dsDay$Rg, VPD = rep(NA_real_, nRow), Temp = dsDay$Temp, VPD0 = 0, fixVPD = TRUE)$GPP
  expect_equal(nrow(gpp), nRow)
  expect_true(all(gpp[, "alpha"] == gpp[1, "alpha"])) # same value not regarding different k
  expect_true(is.finite(gpp[1, "alpha"])) # finite despite NA VPD
  #
  gpp <- lrcFitter$computeLRCGradient(theta0M, Rg = dsDay$Rg, VPD = rep(NA_real_, nRow), Temp = dsDay$Temp, VPD0 = 0)$GPP
  expect_equal(nrow(gpp), nRow)
  expect_true(all(gpp[-iRowK, "alpha"] == gpp[1, "alpha"])) # same value not regarding different k
  expect_true(is.finite(gpp[1, "alpha"])) # finite despite NA VPD
  expect_true(is.na(gpp[iRowK, "alpha"]))
}

test_that("RectangularLRCFitter$predictLRC vector form", {
  lrcFitter <- RectangularLRCFitter()
  test_LRCGRadient_vectorForm(lrcFitter)
  # lrcFitter <- LogisticSigmoidLRCFitter()
})

test_that("NonrectangularLRCFitter$predictLRC vector form", {
  lrcFitter <- NonrectangularLRCFitter()
  test_LRCGRadient_vectorForm(lrcFitter)
})

test_LRCGradient <- function(lrcFitter) {
  # str(ds)
  ds <- dsNEE
  theta0 <- structure(c(0.05, 27.3333395589509, 0.162207578338878, 2.59392002410639, 185, 1.1), .Names = c("k", "beta", "alpha", "RRef", "E0", "logitconv"))
  res <- lrcFitter$computeLRCGradient(theta0, Rg = ds$Rg_f, VPD = ds$VPD_f, Temp = ds$Temp)
  .numDerivLRC <- function(theta, eps = 0.0001, ..., varName = "NEP") {
    ans <- matrix(NA, nrow = length(list(...)[[1]]), ncol = length(theta), dimnames = list(NULL, names(theta)))
    i <- 1L
    for (i in seq_along(theta)) {
      thetaMinus <- theta
      thetaMinus[i] <- theta[i] - eps
      thetaPlus <- theta
      thetaPlus[i] <- theta[i] + eps
      fMinus <- lrcFitter$predictLRC(thetaMinus, ...)[[varName]]
      fPlus <- lrcFitter$predictLRC(thetaPlus, ...)[[varName]]
      ans[, i] <- derivI <- (fPlus - fMinus) / (2 * eps)
    }
    ans
  }
  varName <- "NEP"
  # varName <- "Reco"
  # varName <- "GPP"
  res20 <- .numDerivLRC(varName = varName, theta = theta0, Rg = ds$Rg_f, VPD = ds$VPD_f, Temp = ds$Temp)
  res2 <- res20[, colnames(res[[varName]])]
  expect_true(all(na.omit(abs(res[[varName]] - res2) < 1e-2)))
  # plot( res[[varName]][,2L] ~ res2[,2L])
  # plot( res$NEP[,4L] ~ res2[,4L])
  # plot( resRectGrad$NEP[,3L] ~ res$NEP[,3L])
}

test_that("RectangularLRCFitter$computeLRCGradient matches numerical estimates", {
  lrcFitter <- RectangularLRCFitter()
  test_LRCGradient(lrcFitter)
  # lrcFitter <- LogisticSigmoidLRCFitter()
  # lrcFitter <- NonrectangularLRCFitter()
})


test_that("RHLightResponseCostC", {
  .tmp.reloadDll <- function() {
    library.dynam.unload("REddyProc", file.path(.libPaths()[1], "REddyProc"))
    installPkg()
    library.dynam("REddyProc", "REddyProc", .libPaths()[1])
  }
  #
  dss <- subset(dsNEE, as.POSIXlt(dsNEE$sDateTime)$mday %in% 1:8)
  dssDay <- subset(dss, isDay == TRUE)
  theta <- c(k = 0, beta = 28.6, alpha = 0.18, RRef = 2.87, E0 = 185)
  flux <- dssDay$NEE_f
  sdFlux <- dssDay$NEE_fsd
  betaPrior <- 26
  sdBetaPrior <- 0.3 * betaPrior / sqrt(length(flux))
  parameterPrior <- c(0, betaPrior, 8, 15, 185)
  sdParameterPrior <- c(NA, sdBetaPrior, NA, NA, NA)
  predR <- lrcFitter$predictLRC(
    theta,
    dssDay$Rg_f, dssDay$VPD_f, dssDay$Temp, 10.0, FALSE
  )
  RSSR <- lrcFitter$computeCost(theta, theta, seq_along(theta), flux, sdFlux, parameterPrior, sdParameterPrior,
    Rg = dssDay$Rg_f, VPD = dssDay$VPD_f, Temp = dssDay$Temp
  )
  LRC_CVersion <- RectangularLRCFitterCVersion()
  tmp <- LRC_CVersion$computeCost(theta, theta, seq_along(theta), flux, sdFlux, parameterPrior, sdParameterPrior,
    Rg = dssDay$Rg_f, VPD = dssDay$VPD_f, Temp = dssDay$Temp
  )
  # (predR$NEP - tmp)
  expect_true(tmp - RSSR < 1e-8)
})

.benchmark_RHLightResponseCostC <- function() {
  # require(rbenchmark)
  tmp <- benchmark(
    lrcFitter$computeCost(
      theta[1:4], theta, 1:4, flux, sdFlux, parameterPrior, sdParameterPrior,
      dssDay$Rg_f, dssDay$VPD_f, dssDay$Temp
    ),
    LRC_CVersion$computeCost(
      theta[1:4], theta, 1:4, flux, sdFlux, parameterPrior, sdParameterPrior,
      dssDay$Rg_f, dssDay$VPD_f, dssDay$Temp
    ),
    replications = 10000
  )
  tmp # speedup of only 2 :(
}

test_that("fitLRC", {
  dss <- subset(dsNEE, as.POSIXlt(dsNEE$sDateTime)$mday %in% 1:8)
  dssDay <- subset(dss, isDay == TRUE)
  dssNight <- subset(dss, isNight == TRUE)
  dsDay <- data.frame(NEE = dssDay$NEE_f, sdNEE = dssDay$NEE_fsd, Rg = dssDay$Rg_f, Temp = dssDay$Temp, VPD = dssDay$VPD_f)
  # lrcFitter <- RectangularLRCFitter()
  res <- resRect <- lrcFitter$fitLRC(dsDay,
    E0 = 185, sdE0 = .05 * 185, RRefNight = mean(dssNight$NEE_f, na.rm = TRUE), lastGoodParameters = NA_real_,
    controlGLPart = partGLControl(nBootUncertainty = 10L)
  )
  parNames <- as.vector(lrcFitter$getParameterNames())
  expect_equal(names(res$thetaOpt), parNames)
  expect_equal(names(res$thetaInitialGuess), parNames)
  expect_equal(colnames(res$covParms), parNames)
  expect_true(!all(res$covParms["E0", 1:4] == 0))
  # testing Lasslop compliency: different priors, covariance from fit
  res <- lrcFitter$fitLRC(dsDay,
    E0 = 185, sdE0 = .05 * 185, RRefNight = mean(dssNight$NEE_f, na.rm = TRUE), lastGoodParameters = NA_real_,
    controlGLPart = partGLControl(nBootUncertainty = 0L, isLasslopPriorsApplied = TRUE)
  )
  expect_true(all(res$covParms["E0", 1:4] == 0)) # without bootstrap no covariance wiht other parameters
  # dput(res$opt.parms.V)
  .tmp.plot <- function() {
    dsDay <- dsDay[order(dsDay$Rg), ]
    plot(-NEE ~ Rg, dsDay) # NEE negative?
    p <- res$thetaOpt
    pred <- lrcFitter$predictLRC(p, Rg = dsDay$Rg, VPD = dsDay$VPD, Temp = dsDay$Temp)
    lines(pred$NEP ~ dsDay$Rg)
  }
  # testing increasing number of bootstrap samples
  .tmp.f <- function() {
    (res60 <- lrcFitter$fitLRC(dsDay,
      E0 = 185, sdE0 = .05 * 185, RRefNight = mean(dssNight$NEE_f, na.rm = TRUE),
      controlGLPart = partGLControl(nBootUncertainty = 100L)
    ))
  }
})

#----------------------------- Logistic Sigmoid
test_that("LogisticSigmoidLRCFitter$computeLRCGradient matches numerical estimates", {
  # lrcFitter <- RectangularLRCFitter()
  lrcFitter <- LogisticSigmoidLRCFitter()
  # lrcFitter <- NonrectangularLRCFitter()
  test_LRCGradient(lrcFitter)
})

test_that("fitLRC_LogisticSigmoid", {
  dss <- subset(dsNEE, as.POSIXlt(dsNEE$sDateTime)$mday %in% 1:8)
  dssDay <- subset(dss, isDay == TRUE)
  dssNight <- subset(dss, isNight == TRUE)
  dsDay <- data.frame(NEE = dssDay$NEE_f, sdNEE = dssDay$NEE_fsd, Rg = dssDay$Rg_f, Temp = dssDay$Temp, VPD = dssDay$VPD_f)
  LRCl <- LogisticSigmoidLRCFitter()
  res <- resSig <- LRCl$fitLRC(dsDay,
    E0 = 185, sdE0 = .05 * 185, RRefNight = mean(dssNight$NEE_f, na.rm = TRUE), lastGoodParameters = NA_real_,
    controlGLPart = partGLControl(nBootUncertainty = 10L)
  )
  expect_true(!all(res$covParms["E0", 1:4] == 0))
  # testing Lasslop compliency: different priors, covariance from fit
  res <- LRCl$fitLRC(dsDay,
    E0 = 185, sdE0 = .05 * 185, RRefNight = mean(dssNight$NEE_f, na.rm = TRUE), lastGoodParameters = NA_real_,
    controlGLPart = partGLControl(nBootUncertainty = 0L, isLasslopPriorsApplied = TRUE)
  )
  expect_true(all(res$covParms["E0", 1:4] == 0)) # without bootstrap no covariance wiht other parameters
  # dput(res$opt.parms.V)
  .tmp.plot <- function() {
    # dsDay <- list(...)$dsDay
    dsDay <- dsDay[order(dsDay$Rg), ]
    plot(-NEE ~ Rg, dsDay) # NEE negative?
    p <- res$thetaOpt
    # p <- resOpt$theta
    # LRCn$trace(predictLRC, recover); #LRCn$trace(predictLRC)
    pred <- LRCl$predictLRC(p, Rg = dsDay$Rg, VPD = dsDay$VPD, Temp = dsDay$Temp)
    lines(pred$NEP ~ dsDay$Rg)
    plot(I(pred$NEP + dsDay$NEE) ~ dsDay$Rg)
    abline(0, 0) # inspect residuals
  }
  # testing increasing number of bootstrap samples
  .tmp.f <- function() {
    (res60 <- LRCl$fitLRC(dsDay,
      E0 = 185, sdE0 = .05 * 185, RRefNight = mean(dssNight$NEE_f, na.rm = TRUE), lastGoodParameters = NA_real_,
      controlGLPart = partGLControl(nBootUncertainty = 100L)
    ))
  }
})


#----------------------------- Nonrectangular
test_that("NonrectangularLRCFitter$computeGPPGradient matches numerical estimates", {
  # str(ds)
  ds <- dsNEE
  theta0 <- structure(c(0, 27.3333395589509, 0.162207578338878, 2.59392002410639, 185, 1.1), .Names = c("k", "beta", "alpha", "RRef", "E0", "logitconv"))
  lrcFitter <- NonrectangularLRCFitter()
  res <- lrcFitter$computeGPPGradient(Rg = ds$Rg_f, Amax = theta0["beta"], alpha = theta0["alpha"], logitconv = theta0["logitconv"])
  .numDerivGPP <- function(theta, eps = 0.00001, Rg, ...) {
    GPPParNames <- c("beta", "alpha", "logitconv")
    ans <- matrix(NA, nrow = length(list(...)[[1]]), ncol = length(GPPParNames), dimnames = list(NULL, GPPParNames))
    i <- "beta"
    for (i in GPPParNames) {
      thetaMinus <- theta
      thetaMinus[i] <- theta[i] - eps
      thetaPlus <- theta
      thetaPlus[i] <- theta[i] + eps
      fMinus <- lrcFitter$predictGPP(Rg = Rg, Amax = thetaMinus["beta"], alpha = thetaMinus["alpha"], conv = invlogit(thetaMinus["logitconv"]))
      fPlus <- lrcFitter$predictGPP(Rg = Rg, Amax = thetaPlus["beta"], alpha = thetaPlus["alpha"], conv = invlogit(thetaPlus["logitconv"]))
      ans[, i] <- derivI <- (fPlus - fMinus) / (2 * eps)
    }
    ans
  }
  res2 <- .numDerivGPP(theta = theta0, Rg = ds$Rg_f, VPD = ds$VPD_f, Temp = ds$Temp)
  expect_true(all(na.omit(abs(res - res2) < 1e-2)))
  # plot( res[,2L] ~ res2[,2L])
})

test_that("NonrectangularLRCFitter$computeLRCGradient matches numerical estimates", {
  # lrcFitter <- RectangularLRCFitter()
  # lrcFitter <- LogisticSigmoidLRCFitter()
  lrcFitter <- NonrectangularLRCFitter()
  test_LRCGradient(lrcFitter)
})


test_that("fitLRC_Nonrectangular", {
  dss <- subset(dsNEE, as.POSIXlt(dsNEE$sDateTime)$mday %in% 1:8)
  dssDay <- subset(dss, isDay == TRUE)
  dssNight <- subset(dss, isNight == TRUE)
  dsDay <- data.frame(NEE = dssDay$NEE_f, sdNEE = dssDay$NEE_fsd, Rg = dssDay$Rg_f, Temp = dssDay$Temp, VPD = dssDay$VPD_f)
  # loadPkg()
  LRCn <- NonrectangularLRCFitter()
  res <- resNonrect <- LRCn$fitLRC(dsDay,
    E0 = 185, sdE0 = .05 * 185, RRefNight = mean(dssNight$NEE_f, na.rm = TRUE), lastGoodParameters = NA_real_,
    controlGLPart = partGLControl(nBootUncertainty = 10L)
  )
  expect_true(!all(res$covParms["E0", 1:4] == 0))
  # testing Lasslop compliency: different priors, covariance from fit
  res <- LRCn$fitLRC(dsDay,
    E0 = 185, sdE0 = .05 * 185, RRefNight = mean(dssNight$NEE_f, na.rm = TRUE), lastGoodParameters = NA_real_,
    controlGLPart = partGLControl(nBootUncertainty = 0L, isLasslopPriorsApplied = TRUE)
  )
  expect_true(all(res$covParms["E0", 1:4] == 0)) # without bootstrap no covariance wiht other parameters
  # dput(res$opt.parms.V)
  .tmp.plot <- function() {
    # dsDay <- list(...)$dsDay
    dsDay <- dsDay[order(dsDay$Rg), ]
    plot(-NEE ~ Rg, dsDay) # NEE negative?
    p <- res$thetaOpt
    # p <- resOpt$theta
    # p <- theta0Adj
    # p <- theta0
    # LRCn$trace(predictLRC, recover); #LRCn$trace(predictLRC)
    pred <- LRCn$predictLRC(p, Rg = dsDay$Rg, VPD = dsDay$VPD, Temp = dsDay$Temp)
    lines(pred$NEP ~ dsDay$Rg)
    plot(I(pred$NEP + dsDay$NEE) ~ dsDay$Rg)
    abline(0, 0) # inspect residuals
  }
  # testing increasing number of bootstrap samples
  .tmp.f <- function() {
    (res60 <- LRCn$fitLRC(dsDay,
      E0 = 185, sdE0 = .05 * 185, RRefNight = mean(dssNight$NEE_f, na.rm = TRUE),
      controlGLPart = partGLControl(nBootUncertainty = 100L)
    ))
    sqrt(diag(res60$covParms)) / res60$thetaOpt
  }
})

### ------------------------------ error on inverting hessian ------
args <- list(
  dsDay = .dsDay,
  E0 = 255.598188786627,
  sdE0 = 63.36264,
  RRefNight = 3.47679878519779,
  controlGLPart = list(
    LRCFitConvergenceTolerance = 0.001, nLRCFitConvergenceTolerance = 0.001,
    nBootUncertainty = 0L, minNRecInDayWindow = 10L, isAssociateParmsToMeanOfValids = TRUE,
    isLasslopPriorsApplied = TRUE, isUsingLasslopQualityConstraints = FALSE,
    isSdPredComputed = FALSE, isFilterMeteoQualityFlag = FALSE,
    isBoundLowerNEEUncertainty = TRUE, fixedTRefAtNightTime = NA,
    isExtendTRefWindow = TRUE, smoothTempSensEstimateAcrossTime = TRUE,
    isNeglectPotRadForNight = FALSE, isNeglectVPDEffect = FALSE,
    isRefitMissingVPDWithNeglectVPDEffect = FALSE, fixedTempSens = structure(
      list(
        E0 = NA_real_, sdE0 = NA_real_, RRef = NA_real_
      ),
      class = "data.frame", row.names = c(NA, -1L)
    ),
    replaceMissingSdNEEParms = c(perc = 0.2, minSd = 0.7),
    neglectNEEUncertaintyOnMissing = FALSE, minPropSaturation = NA,
    useNightimeBasalRespiration = FALSE
  ),
  lastGoodParameters = c(
    k = 0.939127165902557, beta = 3.71354808845461, alpha = 0.0382004665603087,
    RRef = 1.05034157955867, E0 = 41.1118669819426
  )
)

test_that("error on inverting hessian", {
  lrcFitter <- RectangularLRCFitter$new()
  # lrcFitter$trace(fitLRC, browser)
  # lrcFitter$trace(optimLRCBounds, browser)
  ans <- do.call(lrcFitter$fitLRC, args)
  expect_true(all(is.na(ans$thetaOpt)))
  # expect_equal(ans$convergence, 1006) # now with this data isFixedVPD is used
})
