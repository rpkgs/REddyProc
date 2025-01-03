
<!-- 
README.md is generated from README.Rmd. Please edit that file
#knitr::knit("README.Rmd") 
rmarkdown::render("README.Rmd") 
maybe clear cache before
-->
<!-- badges: start -->

<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/REddyProc)](http://cran.r-project.org/package=REddyProc) -->
[![R-CMD-check](https://github.com/EarthyScience/REddyProc/workflows/R-CMD-check/badge.svg)](https://github.com/EarthyScience/REddyProc/actions)
[![codecov](https://codecov.io/gh/rpkgs/REddyProc/branch/master/graph/badge.svg)](https://app.codecov.io/gh/rpkgs/REddyProc/tree/master/R)
<!-- badges: end -->

## Overview

`REddyProc` package supports processing (half)hourly data from
Eddy-Covariance sensors.

There is an online-formular to use the functionality of the package
including description at
<https://www.bgc-jena.mpg.de/bgi/index.php/Services/REddyProcWeb>.

## Installation

``` r
# Release stable version from CRAN
install.packages("REddyProc")

# The development version from GitHub using devtools:
# install.packages("devtools")
devtools::install_github("EarthyScience/REddyProc")
```

The REddyProc~package requires a quite recent versions of the tidyverse
packages. On encountering problems on installations with older versions
should run the following code before installing REddyProc.

``` r
install.packages("tidyverse")
update.packages(oldPkgs="dplyr")
```

## Usage

A simple example performs Lookuptable-based gapfilling of
Net-Ecosystem-Exchange (NEE) and plotting a fingerprint plot of the
filled values.

``` r
library(REddyProc)
#+++ Input data from csv (example needs to be downloaded)
examplePath <- getExamplePath('Example_DETha98.txt', isTryDownload = TRUE)
if (length(examplePath)) {
  EddyData <- fLoadTXTIntoDataframe(examplePath)
} else {
  warning(
      "Could not find example text data file."
      ," In order to execute this example code,"
      ," please, allow downloading it from github. " 
      ," Type '?getExamplePath' for more information.")
  # using RData version distributed with the package instead
  EddyData <- Example_DETha98
}
#+++ If not provided, calculate VPD from Tair and rH
EddyData$VPD <- fCalcVPDfromRHandTair(EddyData$rH, EddyData$Tair)
#+++ Add time stamp in POSIX time format
EddyDataWithPosix <- EddyData %>% 
  filterLongRuns("NEE") %>% 
  fConvertTimeToPosix('YDH', Year = 'Year', Day = 'DoY', Hour = 'Hour')
#+++ Initalize R5 reference class sEddyProc for processing of eddy data
#+++ with all variables needed for processing later
EProc <- sEddyProc$new(
  'DE-Tha', EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar'))
#Location of DE-Tharandt
EProc$sSetLocationInfo(LatDeg = 51.0, LongDeg = 13.6, TimeZoneHour = 1)  
#
#++ Fill NEE gaps with MDS gap filling algorithm (without prior ustar filtering)
EProc$sMDSGapFill('NEE', FillAll = FALSE)
#
#++ Export gap filled and partitioned data to standard data frame
FilledEddyData <- EProc$sExportResults()
#
#++ Example plots of filled data to screen or to directory \plots
EProc$sPlotFingerprintY('NEE_f', Year = 1998)
```

![](README-example-1.png)<!-- -->

Further examples are in
[vignette(useCase)](https://github.com/EarthyScience/REddyProc/blob/master/vignettes/useCase.md)
and
[vignette(DEGebExample)](https://github.com/EarthyScience/REddyProc/blob/master/vignettes/DEGebExample.md)
and further md-files of the [vignettes
directory](https://github.com/EarthyScience/REddyProc/blob/master/vignettes).

<!-- ## Docker images

Docker images are provided that comprise rstudio, rocker/tidyverse, and
REddyProc. There are different version for the latest push to github,
for the version on CRAN and for specific tags starting from 1.1.4.

- EarthyScience/reddyproc:latest  
- EarthyScience/reddyproc_cran
- EarthyScience/reddyproc:`tag`

They are usually run with installed docker by typing at a shell:

    docker run --rm -p 8787:8787 -e PASSWORD=REddyProc <imagename>

Then the loading url `localhost:8787` in a browser window should bring
up RStudio  
(default username is rstudio and password was set to REddyProc), where
you can type the above usage example.

For processing your own files in docker you need to mount local
directories with the [–mount
option](https://docs.docker.com/storage/bind-mounts/), e.g.
`--mount type=bind,source=/home/twutz/devR,target=/home/rstudio/devR -e USERID=$UID` -->

## Reference

The methodology and benchmark of `REddyProc` 1.1.3 is described in the
following paper:

Wutzler, T., Lucas-Moffat, A., Migliavacca, M., Knauer, J., Sickel, K.,
Šigut, L., Menzer, O., and Reichstein, M. (2018): Basic and extensible
post-processing of eddy covariance flux data with REddyProc,
Biogeosciences, 15, 5015-5030,
<https://doi.org/10.5194/bg-15-5015-2018>.
