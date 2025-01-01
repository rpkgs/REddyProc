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
  dsDay = structure(
    list(
      sDateTime = structure(
        c(
          987401700, 987403500,
          987405300, 987407100, 987408900, 987410700, 987412500, 987414300,
          987416100, 987417900, 987419700, 987421500, 987423300, 987425100,
          987426900, 987428700, 987430500, 987432300, 987434100, 987435900,
          987437700, 987488100, 987489900, 987491700, 987493500, 987495300,
          987497100, 987498900, 987511500, 987513300, 987515100, 987516900,
          987518700, 987520500, 987522300, 987524100, 987572700, 987574500,
          987576300, 987578100, 987579900, 987581700, 987583500, 987585300,
          987587100, 987588900, 987590700, 987592500, 987594300, 987596100,
          987597900, 987599700, 987601500, 987603300, 987605100, 987606900,
          987608700, 987610500, 987659100, 987662700, 987664500, 987666300,
          987668100, 987669900, 987671700, 987673500, 987675300, 987677100,
          987678900, 987680700, 987682500, 987684300, 987686100, 987687900,
          987689700, 987691500, 987693300, 987695100, 987696900
        ),
        class = c("POSIXct", "POSIXt"), tzone = "UTC"
      ),
      NEE = c(
        0.35400003194809, 0.532000005245209,
        0.452000021934509, 0.887999951839447, 0.284999996423721, 0.696999967098236,
        0.50900000333786, 0.435000032186508, 0.32600000500679, 0.445999979972839,
        0.397000014781952, 1.02100002765656, 0.50299996137619, 0.234999999403954,
        0.340000003576279, 0.472000002861023, 0.395000010728836, 0.155000001192093,
        0.21299996972084, 0.679999947547913, 0.42399999499321, 0.471999943256378,
        0.187000036239624, -0.318000018596649, 0.64599996805191, 0.204000011086464,
        0.367000013589859, -1.28100001811981, 0.0670000314712524, 1.12699997425079,
        0.709000051021576, 1.16599988937378, 3.26600003242493, 0.282000005245209,
        1.04100000858307, 1.20700001716614, 0.533999979496002, -0.175999999046326,
        -0.256999999284744, -1.32000005245209, 0.188999950885773, 1.66199994087219,
        0.196999996900558, 0.566999971866608, 0.039000004529953, 0.725000023841858,
        -0.650000035762787, -0.583000004291534, 0.441000014543533, 1.19400000572205,
        1.04999995231628, 0.671000003814697, 0.187000006437302, -0.412999987602234,
        -0.306999981403351, 0.785999953746796, 0.638999998569489, 0.738000094890594,
        -0.340000003576279, 1.52699995040894, 0.712999999523163, 0.130999982357025,
        -0.101999998092651, -0.199999988079071, 0.453000009059906, -0.0219999849796295,
        1.14800000190735, -0.293000012636185, 0.880000054836273, 0.678999960422516,
        0.795000016689301, -0.0259999930858612, 0.787999987602234, 1.27799999713898,
        -0.248999997973442, -0.462000012397766, 0.874999940395355, 0.51800000667572,
        0.36300003528595
      ),
      sdNEE = c(
        0.809749957666536, 0.628224972940061,
        0.659803884323685, 0.662579436599356, 0.650053791385915, 0.660612423407973,
        0.610460777640826, 0.542368468171728, 0.594734288978246, 0.579223376405118,
        0.630268911933181, 0.607940657100616, 0.659803884323685, 0.575262505769259,
        0.475729669132864, 0.648430658832156, 0.682795725496793, 0.577201698568239,
        0.698479742960692, 0.71378981486918, 0.606161513437173, 0.786607636006147,
        0.806993072704408, 0.684570302164522, 0.723004813652185, 0.623127617572988,
        0.718289427337163, 0.732915617698706, 0.718289427337163, 0.58019133449497,
        0.685697693495996, 0.668489300649783, 1.07041296286867, 0.398845289016157,
        0.470392568692144, 0.697977377661453, 0.787856461757981, 1.03571321161812,
        0.938838173853171, 0.697215779779815, 0.453336654077684, 0.882031438984157,
        1.16052219915433, 0.555610054123532, 0.401003326590866, 0.578448314947156,
        0.597527646859134, 0.599727434290068, 0.579782790464044, 0.5773945573756,
        0.599727434290068, 0.536770541502652, 0.573246404355542, 1.08054803991486,
        0.733622677330087, 0.655104966509601, 0.752326787415782, 0.529818303846866,
        0.806993072704408, 0.806695090744094, 0.737711291476921, 0.576470661150769,
        0.72864610563661, 0.677348267331516, 0.653979307248866, 0.48473545601523,
        1.01680853513454, 0.981286911668651, 0.466287618250474, 0.331923568560297,
        1.04350607959198, 0.48138096273862, 0.563125213132485, 1.01610379909119,
        0.727864584639489, 0.604706926370824, 0.618848471216863, 0.691699568925105,
        0.461768677852474
      ),
      Temp = c(
        1.20500004291534, 1.60399997234344,
        2.35199999809265, 1.41700005531311, 1.60300004482269, 1.82299995422363,
        1.86800003051758, 2.2409999370575, 2.17499995231628, 2.28500008583069,
        2.125, 2.25399994850159, 2.41599988937378, 2.33800005912781,
        2.21099996566772, 1.99300003051758, 2.09400010108948, 1.99899995326996,
        2.05100011825562, 2.09999990463257, 2.17000007629395, 2.04200005531311,
        1.39800000190735, 1.40199995040894, 1.68700003623962, 1.88300001621246,
        1.87999999523163, 1.52900004386902, 1.86699998378754, 2.59100008010864,
        2.64100003242493, 2.75200009346008, 4.47800016403198, 4.32499980926514,
        4.64300012588501, 3.9210000038147, -1.807000041008, -1.78199994564056,
        -0.637000024318695, 1.91299998760223, 3.66100001335144, 4.02299976348877,
        4.33199977874756, 4.75299978256226, 5.06799983978271, 4.91200017929077,
        6.28700017929077, 6.38600015640259, 6.65299987792969, 6.44199991226196,
        6.57200002670288, 6.60300016403198, 7.10500001907349, 6.75799989700317,
        7.23000001907349, 7.92500019073486, 7.97499990463257, 6.31699991226196,
        1.37300002574921, 1.57200002670288, 1.83299994468689, 2.93600010871887,
        2.7590000629425, 2.99399995803833, 4.08699989318848, 5.32800006866455,
        5.83099985122681, 5.95300006866455, 5.80900001525879, 5.70599985122681,
        5.93200016021729, 6.03800010681152, 7.23099994659424, 6.81300020217896,
        7.46600008010864, 7.76999998092651, 8.625, 8.92500019073486,
        9.10000038146973
      ),
      VPD = c(
        0.200011223554611, 0.205833956599236,
        0.231631278991699, 0.223395854234695, 0.219540312886238, 0.216062813997269,
        0.202776655554771, 0.19389633834362, 0.178689301013947, 0.165691033005714,
        0.15668548643589, 0.158136337995529, 0.15270359814167, 0.144625097513199,
        0.13615371286869, 0.12699282169342, 0.12080705165863, 0.112930752336979,
        0.106267288327217, 0.0995310693979263, 0.0928854197263718, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0252391025424004, 0.307956606149673,
        1.0468270778656, 1.16496980190277, 0.834864974021912, 1.00261807441711,
        1.2131587266922, 1.90804219245911, 2.35116386413574, 2.20005965232849,
        2.20672249794006, 2.39285826683044, 2.51623201370239, 2.50632429122925,
        3.08190584182739, 3.11264896392822, 3.08240079879761, 2.98003649711609,
        3.3766405582428, 4.0762357711792, 4.41115427017212, 4.14953422546387,
        4.072425365448, 4.51607370376587, 4.21010255813599, 3.14566493034363,
        0.573594272136688, 0.568167388439178, 0.453360319137573, 0.430106312036514,
        0.476890712976456, 0.545533657073975, 0.818514168262482, 1.62483441829681,
        1.95072364807129, 2.30293464660645, 2.07698583602905, 2.13552165031433,
        2.14131665229797, 2.20398330688477, 2.91198372840881, 2.90875840187073,
        3.2385835647583, 3.56008434295654, 4.03070735931396, 4.26191759109497,
        4.86756706237793
      ),
      Rg = c(
        19.71875, 105.892852783203, 131.607131958008,
        37.4553565979004, 59.8214302062988, 53.035717010498, 88.5714263916016,
        141.741073608398, 131.473220825195, 157.857147216797, 75.2678604125977,
        89.7767868041992, 129.866073608398, 115.044647216797, 115.044647216797,
        58.0803604125977, 47.7232131958008, 21.1339282989502, 38.4375,
        36.9196395874023, 25.9285717010498, 9.4776782989502, 10.6741065979004,
        39.8660736083984, 83.5714263916016, 111.116065979004, 101.071426391602,
        91.6071395874023, 99.6428527832031, 170.982147216797, 185.178558349609,
        215.982131958008, 418.75, 299.821411132812, 216.607147216797,
        190.267868041992, 14.8169631958008, 37.3660697937012, 110.223213195801,
        225.178573608398, 313.392852783203, 365.178558349609, 436.160705566406,
        508.482147216797, 532.589294433594, 507.589294433594, 588.392883300781,
        655.803588867188, 626.785705566406, 622.321411132812, 654.464294433594,
        644.642883300781, 646.875, 450.892852783203, 477.678558349609,
        463.839294433594, 408.035705566406, 252.991073608398, 14.5982141494751,
        11.9821424484253, 96.7410659790039, 173.4375, 134.866073608398,
        145, 239.732147216797, 318.75, 402.232147216797, 403.571441650391,
        324.107147216797, 303.258911132812, 379.017852783203, 353.571441650391,
        575.446411132812, 376.339294433594, 445.089294433594, 496.428558349609,
        463.839294433594, 402.678558349609, 337.5
      ),
      isDay = c(
        TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      isNight = c(
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
      )
    ),
    row.names = c(
      92749L, 92750L, 92751L, 92752L, 92753L, 92754L,
      92755L, 92756L, 92757L, 92758L, 92759L, 92760L, 92761L, 92762L,
      92763L, 92764L, 92765L, 92766L, 92767L, 92768L, 92769L, 92797L,
      92798L, 92799L, 92800L, 92801L, 92802L, 92803L, 92810L, 92811L,
      92812L, 92813L, 92814L, 92815L, 92816L, 92817L, 92844L, 92845L,
      92846L, 92847L, 92848L, 92849L, 92850L, 92851L, 92852L, 92853L,
      92854L, 92855L, 92856L, 92857L, 92858L, 92859L, 92860L, 92861L,
      92862L, 92863L, 92864L, 92865L, 92892L, 92894L, 92895L, 92896L,
      92897L, 92898L, 92899L, 92900L, 92901L, 92902L, 92903L, 92904L,
      92905L, 92906L, 92907L, 92908L, 92909L, 92910L, 92911L, 92912L,
      92913L
    ), class = "data.frame"
  ),
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
