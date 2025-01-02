ds <- partGLExtractStandardData(dsNEE)
dsTempSens <- partGL_FitNight_E0_RRef(ds, nRecInDay = 48L, controlGLPart = partGLControl())

# lrcFitter <- NonrectangularLRCFitter()
lrcFitter <- RectangularLRCFitter()
resFits <- partGLFitLRCWindows(ds, nRecInDay = 48L, 
  dsTempSens = dsTempSens, lrcFitter = lrcFitter,
  controlGLPart = partGLControl(nBootUncertainty = 10L))
resFits
