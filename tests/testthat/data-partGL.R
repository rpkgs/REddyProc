## VARIABLES: resLRCEx1, resLRCEx1Nonrectangular

# regression result from dput in
resLRCEx1 <- structure(
  list(
    iWindow = 1:4, dayStart = c(1L, 3L, 5L, 7L), dayEnd = c(4L, 6L, 8L, 9L),
    iRecStart = c(1L, 97L, 193L, 289L), iRecEnd = c(192L, 288L, 384L, 428L),
    iCentralRec = c(97L, 193L, 289L, 385L), nValidRec = c(98L, 98L, 75L, 56L),
    iMeanRec = c(97L, 194L, 283L, 385L),
    convergence = c(0L, 0L, 0L, 0L), parms_out_range = c(0L, 0L, 0L, 0L),
    k = c(
      0.15556543911127, 0.100336071150363, 0.0292272454449355, 0.110909241852741
    ),
    beta = c(
      34.2363908612952, 28.7053014210465, 11.965007963525, 35.6061046130658
    ),
    alpha = c(
      0.0408287988682114, 0.0317164138158062, 0.12184364495971, 0.0565989757582912
    ),
    RRef = c(
      1.30588274435898, 2.30342378833466, 4.39970416578445, 3.02347421699605
    ),
    E0 = structure(
      c(62.633343296311, 62.633343296311, 62.633343296311, 62.633343296311),
      .Dim = c(4L, 1L), .Dimnames = list(NULL, NULL)
    ),
    k_sd = c(
      0.0162347703672834, 0.00819228863397153, 0.0036484758449877, 0.0143126637303779
    ),
    beta_sd = c(
      2.6633120785356, 2.53043346467805, 0.480085981083196, 3.86498333963868
    ),
    alpha_sd = c(
      0.00282359981111383, 0.00245807510249875, 0.0189581631026536, 0.00624526974304941
    ),
    RRef_sd = c(
      0.243651749533536, 0.200219643557999, 0.379555581976541, 0.408032180867498
    ),
    E0_sd = structure(
      c(4.25106879296635, 4.22511397480583, 4.22511362599932, 4.46593783063789),
      .Dim = c(4L, 1L), .Dimnames = list(NULL, NULL)
    ),
    GPP2000 = c(
      24.1225750077742, 19.7622684102263, 11.4050231061514, 27.0862112802389
    ),
    isValid = c(TRUE, TRUE, TRUE, TRUE),
    resOpt = list(
      structure(
        list(
          thetaOpt = structure(
            c(
              0.15556543911127, 34.2363908612952, 0.0408287988682114,
              1.30588274435898, 62.633343296311
            ),
            .Names = c("k", "beta", "alpha", "RRef", "E0")
          ),
          iOpt = 1:5, thetaInitialGuess = structure(
            c(
              0.05, 24.3542, 0.1, 5.17329116845493, 62.633343296311
            ),
            .Names = c("k", "beta", "alpha", "RRef", "E0")
          ),
          covParms = structure(
            c(
              0.000263567768878423, 0.0299004485274175, -2.01722085100335e-05,
              -0.00170017321010028, 0, 0.0299004485274174, 7.09323122767363,
              -0.00442973195793744, -0.184157427431281, 0, -2.01722085100335e-05,
              -0.00442973195793744, 7.97271589332207e-06, 0.000574595628371077, 0,
              -0.00170017321010028, -0.18415742743128, 0.000574595628371077,
              0.0593661750507529, 0, 0, 0, 0, 0, 18.0715858825324
            ),
            .Dim = c(5L, 5L), .Dimnames = list(
              structure(c("k", "beta", "alpha", "RRef", "E0"),
                .Names = c("k", "beta", "alpha", "RRef", "E0")
              ),
              structure(c("k", "beta", "alpha", "RRef", "E0"),
                .Names = c("k", "beta", "alpha", "RRef", "E0")
              )
            )
          ),
          convergence = 0L
        ),
        .Names = c(
          "thetaOpt", "iOpt", "thetaInitialGuess", "covParms", "convergence"
        )
      ),
      structure(
        list(
          thetaOpt = structure(
            c(
              0.100336071150363, 28.7053014210465, 0.0317164138158062,
              2.30342378833466, 62.633343296311
            ),
            .Names = c("k", "beta", "alpha", "RRef", "E0")
          ),
          iOpt = 1:5,
          thetaInitialGuess = structure(
            c(0.05, 24.8403, 0.1, 5.14770049230269, 62.633343296311),
            .Names = c("k", "beta", "alpha", "RRef", "E0")
          ),
          covParms = structure(
            c(
              6.71135930622992e-05, 0.0178801034106948, -1.74280275277554e-05,
              -0.00134394852514135, 0, 0.0178801034106948, 6.40309351916256,
              -0.00493560933704443, -0.271156219332448, 0, -1.74280275277554e-05,
              -0.00493560933704443, 6.04213320952425e-06, 0.000418391213333105,
              0, -0.00134394852514135, -0.271156219332449, 0.000418391213333105,
              0.040087905666492, 0, 0, 0, 0, 0, 17.8515881000995
            ),
            .Dim = c(5L, 5L), .Dimnames = list(
              structure(c("k", "beta", "alpha", "RRef", "E0"),
                .Names = c("k", "beta", "alpha", "RRef", "E0")
              ),
              structure(c("k", "beta", "alpha", "RRef", "E0"),
                .Names = c("k", "beta", "alpha", "RRef", "E0")
              )
            )
          ),
          convergence = 0L
        ),
        .Names = c(
          "thetaOpt", "iOpt", "thetaInitialGuess", "covParms", "convergence"
        )
      ),
      structure(
        list(
          thetaOpt = structure(
            c(
              0.0292272454449355, 11.965007963525, 0.12184364495971,
              4.39970416578445, 62.633343296311
            ),
            .Names = c("k", "beta", "alpha", "RRef", "E0")
          ),
          iOpt = 1:5,
          thetaInitialGuess = structure(
            c(0.05, 24.1194, 0.1, 5.14770049230269, 62.633343296311),
            .Names = c("k", "beta", "alpha", "RRef", "E0")
          ),
          covParms = structure(
            c(
              1.33113759914587e-05, 0.000439566354933103, -5.94748188903169e-05,
              -0.00111894361237909, 0, 0.000439566354933103, 0.230482549232615,
              -0.000116871653668728, 0.053260722726946, 0, -5.94748188903169e-05,
              -0.000116871653668727, 0.000359411948226816, 0.00657062463616557,
              0, -0.00111894361237909, 0.0532607227269459, 0.00657062463616558,
              0.144062439809551, 0, 0, 0, 0, 0, 17.8515851526051
            ),
            .Dim = c(5L, 5L), .Dimnames = list(
              structure(c("k", "beta", "alpha", "RRef", "E0"),
                .Names = c("k", "beta", "alpha", "RRef", "E0")
              ),
              structure(c("k", "beta", "alpha", "RRef", "E0"),
                .Names = c("k", "beta", "alpha", "RRef", "E0")
              )
            )
          ),
          convergence = 0L
        ),
        .Names = c(
          "thetaOpt", "iOpt", "thetaInitialGuess", "covParms", "convergence"
        )
      ),
      structure(
        list(
          thetaOpt = structure(
            c(
              0.110909241852741, 35.6061046130658, 0.0565989757582912,
              3.02347421699605, 62.633343296311
            ),
            .Names = c("k", "beta", "alpha", "RRef", "E0")
          ),
          iOpt = 1:5,
          thetaInitialGuess = structure(
            c(0.05, 21.665, 0.1, 5.21501814931109, 62.633343296311),
            .Names = c("k", "beta", "alpha", "RRef", "E0")
          ),
          covParms = structure(
            c(
              0.000204852343058874, 0.0495123146513288, -6.84478020597351e-05,
              -0.00386015942484132, 0, 0.0495123146513288, 14.9380962156846,
              -0.0178713249608098, -0.786540464781837, 0, -6.8447802059735e-05,
              -0.0178713249608097, 3.90033941634485e-05, 0.00226964418347153,
              0, -0.00386015942484131, -0.786540464781834, 0.00226964418347153,
              0.166490260623486, 0, 0, 0, 0, 0, 19.9446007071226
            ),
            .Dim = c(5L, 5L), .Dimnames = list(
              structure(
                c("k", "beta", "alpha", "RRef", "E0"),
                .Names = c("k", "beta", "alpha", "RRef", "E0")
              ),
              structure(
                c("k", "beta", "alpha", "RRef", "E0"),
                .Names = c("k", "beta", "alpha", "RRef", "E0")
              )
            )
          ),
          convergence = 0L
        ),
        .Names = c(
          "thetaOpt", "iOpt", "thetaInitialGuess", "covParms", "convergence"
        )
      )
    ),
    E0_bootstrap_sd = c(
      4.25106879296635, 4.22511397480583, 4.22511362599932, 4.46593783063789
    ),
    RRef_night = c(
      5.17329116845493, 5.14770049230269, 5.14770049230269, 5.21501814931109
    )
  ),
  .Names = c(
    "iWindow", "dayStart", "dayEnd", "iRecStart", "iRecEnd", "iCentralRec",
    "nValidRec", "iMeanRec", "convergence", "parms_out_range", "k", "beta",
    "alpha", "RRef", "E0", "k_sd", "beta_sd", "alpha_sd", "RRef_sd", "E0_sd",
    "GPP2000", "isValid", "resOpt", "E0_bootstrap_sd", "RRef_night"
  ),
  row.names = c(NA, -4L),
  class = c("tbl_df", "tbl", "data.frame")
)


resLRCEx1Nonrectangular <- structure(
  list(
    iWindow = 1:4, dayStart = c(1L, 3L, 5L, 7L), dayEnd = c(4L, 6L, 8L, 9L),
    iRecStart = c(1L, 97L, 193L, 289L), iRecEnd = c(192L, 288L, 384L, 428L),
    iCentralRec = c(97L, 193L, 289L, 385L), nValidRec = c(98L, 98L, 75L, 56L),
    iMeanRec = c(97L, 194L, 283L, 385L), convergence = c(1005L, 1005L, 1005L, 0L),
    E0 = structure(
      c(
        62.633343296311, 62.633343296311, 62.633343296311, 62.633343296311
      ),
      .Dim = c(4L, 1L), .Dimnames = list(NULL, NULL)
    ),
    E0_sd = structure(
      c(4.25106879296635, 4.22511397480583, 4.22511362599932, 4.46593783063789),
      .Dim = c(4L, 1L), .Dimnames = list(NULL, NULL)
    ),
    RRefNight = c(5.17329116845493, 5.14770049230269, 5.14770049230269, NA),
    parms_out_range = c(NA, NA, NA, 0L), k = c(NA, NA, NA, 0.15301818456226),
    beta = c(NA, NA, NA, 49.3407957035085),
    alpha = c(NA, NA, NA, 0.0444237102808321),
    RRef = c(NA, NA, NA, 2.35368171175251),
    logitconv = c(NA, NA, NA, -3.10163338239329),
    k_sd = c(NA, NA, NA, 0.0197088986436212),
    beta_sd = c(NA, NA, NA, 7.58853578454936),
    alpha_sd = c(NA, NA, NA, 0.00432079403235484),
    RRef_sd = c(NA, NA, NA, 0.350880271587216),
    logitconv_sd = c(NA, NA, NA, 1.92633463362266),
    GPP2000 = c(NA, NA, NA, 32.0432128516724),
    isValid = c(NA, NA, NA, TRUE),
    resOpt = list(NULL, NULL, NULL, structure(
      list(
        thetaOpt = structure(
          c(
            0.15301818456226, 49.3407957035085, 0.0444237102808321, 2.35368171175251,
            62.633343296311, -3.10163338239329
          ),
          .Names = c("k", "beta", "alpha", "RRef", "E0", "logitconv")
        ),
        iOpt = c(1, 2, 3, 4, 6, 5),
        thetaInitialGuess = structure(
          c(
            0.05, 21.665, 0.1, 5.21501814931109, 62.633343296311, 1.09861228866811
          ),
          .Names = c("k", "beta", "alpha", "RRef", "E0", "logitconv")
        ),
        covParms = structure(
          c(
            0.000388440685744535, 0.13746635747882, -6.51906051810543e-05,
            -0.00449519721692684, 0, 0.00207560732637497, 0.13746635747882,
            57.5858753533862, -0.0238832979670161, -1.36443384766771, 0,
            -1.52788892519458, -6.51906051810543e-05, -0.023883297967016,
            1.86692610700332e-05, 0.00131289852541519, 0, -0.00181474371141437,
            -0.00449519721692684, -1.36443384766771, 0.00131289852541519,
            0.123116964989119, 0, -0.0864363703508266, 0, 0, 0, 0, 19.9446007071226,
            0, 0.00207560732637499, -1.52788892519458, -0.00181474371141437,
            -0.0864363703508268, 0, 3.71076512069413
          ),
          .Dim = c(6L, 6L), .Dimnames = list(
            structure(
              c(
                "k", "beta", "alpha", "RRef", "E0", "logitconv"
              ),
              .Names = c("k", "beta", "alpha", "RRef", "E0", "logitconf")
            ),
            structure(
              c(
                "k", "beta", "alpha", "RRef", "E0", "logitconv"
              ),
              .Names = c("k", "beta", "alpha", "RRef", "E0", "logitconf")
            )
          )
        ),
        convergence = 0L
      ),
      .Names = c(
        "thetaOpt", "iOpt", "thetaInitialGuess", "covParms", "convergence"
      )
    )),
    E0_bootstrap_sd = c(
      4.25106879296635, 4.22511397480583, 4.22511362599932, 4.46593783063789
    ),
    RRef_night = c(
      5.17329116845493, 5.14770049230269, 5.14770049230269, 5.21501814931109
    )
  ),
  .Names = c(
    "iWindow", "dayStart", "dayEnd", "iRecStart", "iRecEnd", "iCentralRec",
    "nValidRec", "iMeanRec", "convergence", "E0", "E0_sd", "RRefNight",
    "parms_out_range", "k", "beta", "alpha", "RRef", "logitconv", "k_sd",
    "beta_sd", "alpha_sd", "RRef_sd", "logitconv_sd", "GPP2000", "isValid",
    "resOpt", "E0_bootstrap_sd", "RRef_night"
  ),
  row.names = c(NA, -4L),
  class = c("tbl_df", "tbl", "data.frame")
)
